@echo off
cd /d "%~dp0"
echo.
echo  ========================================
echo   本地测试服务器
echo  ========================================
echo.
echo   密码入口:  http://localhost:8080/
echo   生日站:    http://localhost:8080/7749658f7c83bce1c44ce6f4a514c8c40bb57879/birthday/
echo   入口密码:  082917
echo.
echo   修改文件后刷新浏览器即可，无需 git push
echo   按 Ctrl+C 停止服务器
echo.
python -m http.server 8080
