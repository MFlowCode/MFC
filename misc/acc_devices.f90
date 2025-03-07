program main

    use openacc

    implicit none

    integer :: i
    integer :: num_devices
    integer(acc_device_kind) :: devtype

    character*256 :: devvendor, devname
    real :: devmem
    
    devtype     = acc_get_device_type()
    num_devices = acc_get_num_devices(devtype)

    print '("Found "I0" device(s):")', num_devices
    
    do i = 1, num_devices
        devmem = real(acc_get_property(i, devtype, acc_property_memory)) / 1073741824
        call acc_get_property_string(i, devtype, acc_property_vendor, devvendor)
        call acc_get_property_string(i, devtype, acc_property_name, devname)
        print '(" - "I3" : "A" | "F0.2" GB | "A"")', i, devvendor, devmem, devname
    end do

end program main
