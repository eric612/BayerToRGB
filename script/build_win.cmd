@echo off
@setlocal EnableDelayedExpansion

if NOT EXIST build mkdir build
pushd build

if NOT DEFINED CMAKE_CONFIG set CMAKE_CONFIG=Release

cmake -DCMAKE_GENERATOR="Visual Studio 14 2015 Win64" ^
	  -DCMAKE_BUILD_TYPE:STRING=%CMAKE_CONFIG% ^
	  -DOpenCV_Enable=OFF ^
      ::-DCMAKE_PREFIX_PATH="OpenCV" ^
      "%~dp0\.."
:: Build the library and tools
cmake --build . --config %CMAKE_CONFIG%	 
if ERRORLEVEL 1 (
  echo ERROR: Build failed
  exit /b 1
) 


cd ..	  
popd
@endlocal