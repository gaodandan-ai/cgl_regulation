[Setup]
AppName=Cgl Regulation Explorer
AppVersion=1.3.0
DefaultDirName={autopf}\CglRegulation
DefaultGroupName=CglRegulation
OutputDir=dist
OutputBaseFilename=cgl_setup
Compression=lzma
SolidCompression=yes
DisableProgramGroupPage=no
PrivilegesRequired=lowest
SetupIconFile=icon.ico

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Additional shortcuts:"; Flags: unchecked
Name: "startmenuicon"; Description: "Create a &Start Menu shortcut"; GroupDescription: "Additional shortcuts:"; Flags: checkedonce

[Files]
Source: "dist\cgl_regulation.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "icon.ico"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{userdesktop}\Cgl Regulation Explorer"; Filename: "{app}\cgl_regulation.exe"; WorkingDir: "{app}"; IconFilename: "{app}\icon.ico"; Tasks: desktopicon
Name: "{group}\Cgl Regulation Explorer"; Filename: "{app}\cgl_regulation.exe"; WorkingDir: "{app}"; IconFilename: "{app}\icon.ico"; Tasks: startmenuicon
