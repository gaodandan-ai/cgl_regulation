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

[Files]
Source: "dist\cgl_regulation.exe"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{userdesktop}\C. glutamicum 调控网络浏览器"; Filename: "{app}\cgl_regulation.exe"; WorkingDir: "{app}"
Name: "{group}\C. glutamicum 调控网络浏览器"; Filename: "{app}\cgl_regulation.exe"; WorkingDir: "{app}"
