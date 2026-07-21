[Setup]
AppName=Cgl Regulation Explorer
AppVersion=1.0.0
DefaultDirName={autopf}\CglRegulation
DefaultGroupName=CglRegulation
OutputDir=dist
OutputBaseFilename=cgl_setup
Compression=lzma
SolidCompression=yes
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
SetupIconFile=icon.ico

[Files]
Source: "dist\cgl_regulation.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "icon.ico"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{userdesktop}\Cgl Regulation Explorer"; Filename: "{app}\cgl_regulation.exe"; WorkingDir: "{app}"; IconFilename: "{app}\icon.ico"
Name: "{group}\Cgl Regulation Explorer"; Filename: "{app}\cgl_regulation.exe"; WorkingDir: "{app}"; IconFilename: "{app}\icon.ico"
