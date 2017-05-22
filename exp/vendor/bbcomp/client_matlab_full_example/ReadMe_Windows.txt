If you have some problems loading bbcomp.dll with loadlibrary function under Windows, please check the following.
1. Please make sure you use correct library files, e.g., for 64bit MatLab use 64bit dll files available in bbcomp-windows-64bit.zip. 
2. If you get errors from loadlibrary, please select a compatible compiler by calling "mex -setup" in MatLab.
3. If no compiler is installed (often the case for 64bit MatLab), 
   please see http://mathworks.com/support/compilers/R2012b/win64.html or the same link for your version of MatLab.
   It says that you may/should install Microsoft Windows SDK 7.1 and .NET Framework 4.0. 
   Please CAREFULLY read and FOLLOW instructions available at http://mathworks.com/matlabcentral/answers/95039-why-does-the-sdk-7-1-installation-fail-with-an-installation-failed-message-on-my-windows-system
   SDK 7.1 is available at http://www.microsoft.com/en-us/download/details.aspx?id=8279