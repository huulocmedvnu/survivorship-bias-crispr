"""Seekable read-only file over HTTP range requests, IPv4 only (the HF CDN is unreachable over IPv6 here),
with retries and backoff so a dropped connection does not end the job."""
import io, time
import urllib3.util.connection as _c
_c.HAS_IPV6 = False                      # force IPv4 for requests/urllib3
import requests

S = requests.Session()


class HTTPFile(io.RawIOBase):
    def __init__(self, url, retries=8):
        self.retries = retries
        r = self._req("HEAD", url, allow_redirects=True)
        self.url = r.url                                   # signed CDN URL after the redirect
        self.size = int(r.headers["Content-Length"])
        self.pos = 0
        self.bytes_read = 0

    def _req(self, method, url, **kw):
        for k in range(self.retries):
            try:
                r = S.request(method, url, timeout=60, **kw)
                if r.status_code in (200, 206):
                    return r
            except requests.RequestException:
                pass
            time.sleep(min(60, 2 ** k))
        raise IOError(f"failed after {self.retries} tries: {url}")

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        r = self._req("GET", self.url, headers={"Range": f"bytes={self.pos}-{self.pos + n - 1}"})
        b = r.content
        self.pos += len(b); self.bytes_read += len(b)
        return b

    def readinto(self, buf):
        b = self.read(len(buf)); buf[:len(b)] = b; return len(b)
