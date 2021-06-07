function mountpoint = mount_drive(volume,user,server,mountname)
% mountpoint = mount_drive(volume,user,server,mountname)
% e.g. 
% volume = 'DATA';
% user = 'zeilon';
% server = 'eilon.ldeo.columbia.edu';

if nargin< 4
    mountname=volume;
end

pwd = passcode;

eval(sprintf('! mkdir /Volumes/%s',mountname))
eval(sprintf('! mount_afp "afp://%s:%s@%s/%s" "/Volumes/%s"',user,pwd,server,volume,mountname))

mountpoint = ['/Volumes/',mountname];

end