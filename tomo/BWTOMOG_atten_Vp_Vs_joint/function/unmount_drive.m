function unmount_drive(mountpoint)
% unmount_drive(mountpoint)
eval(sprintf('! umount %s',mountpoint))
end