#:def HardcodedDeallocation()
    if (allocated(stored_values)) then
        @:DEALLOCATE(stored_values)
        @:DEALLOCATE(x_coords)
    end if

    if (allocated(y_coords)) then
        @:DEALLOCATE(y_coords)
    end if

    if (allocated(ih)) then
        @:DEALLOCATE(ih)
    end if
#:enddef
