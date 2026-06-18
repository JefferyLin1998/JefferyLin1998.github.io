const http = require("http");
const fs = require("fs");
const path = require("path");

const root = "D:/JefferyLin1998.github.io";
const types = {
  ".html": "text/html;charset=utf-8",
  ".js": "text/javascript",
  ".css": "text/css",
  ".json": "application/json",
  ".ico": "image/x-icon",
};

http
  .createServer((req, res) => {
    let p = decodeURIComponent(req.url.split("?")[0]);
    let fp = path.join(root, p);

    try {
      const stat = fs.statSync(fp);
      if (stat.isDirectory()) {
        fp = path.join(fp, "index.html");
      }
    } catch (e) {
      // ignore, will 404 below
    }

    fs.readFile(fp, (e, d) => {
      if (e) {
        res.writeHead(404);
        res.end("404 Not Found: " + p);
        return;
      }
      res.writeHead(200, {
        "Content-Type": types[path.extname(fp)] || "application/octet-stream",
      });
      res.end(d);
    });
  })
  .listen(8081, () => console.log("Server running at http://localhost:8081/"));
