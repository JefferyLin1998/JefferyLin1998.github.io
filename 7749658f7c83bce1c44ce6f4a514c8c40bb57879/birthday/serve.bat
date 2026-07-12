@echo off
chcp 65001 >nul
cd /d "%~dp0..\.."

set "PREVIEW_URL=http://localhost:8080/7749658f7c83bce1c44ce6f4a514c8c40bb57879/birthday/"

echo.
echo  ========================================
echo   Birthday 本地预览
echo  ========================================
echo.
echo   生日站直达: %PREVIEW_URL%
echo   密码入口:   http://localhost:8080/
echo   入口密码:   082917
echo.
echo   修改文件后刷新浏览器即可
echo   按 Ctrl+C 停止服务器
echo.

start "" "%PREVIEW_URL%"
python -m http.server 8080
