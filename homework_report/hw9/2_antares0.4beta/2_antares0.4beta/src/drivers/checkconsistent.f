        subroutine checkconsistent()
        implicit none
        include 'fluid.h'
        include 'parameters'
        
        ! check dimension
        if (D1D) then
                if (D2D) then
                        write(*,*)'error: can only be 1D, 2D or 3D '
                        write(*,*)'       checking D1D,D2D,D3D     '
                        stop
                endif
                if (D3D) then
                        write(*,*)'error: can only be 1D, 2D or 3D '
                        write(*,*)'       checking D1D,D2D,D3D     '
                        stop
                endif
        elseif(D2D) then
                if (D1D) then
                        write(*,*)'error: can only be 1D, 2D or 3D '
                        write(*,*)'       checking D1D,D2D,D3D     '
                        stop
                endif
                if (D3D) then
                        write(*,*)'error: can only be 1D, 2D or 3D '
                        write(*,*)'       checking D1D,D2D,D3D     '
                        stop
                endif
        elseif(D3D) then
                if (D1D) then
                        write(*,*)'error: can only be 1D, 2D or 3D '
                        write(*,*)'       checking D1D,D2D,D3D     '
                        stop
                endif
                if (D2D) then
                        write(*,*)'error: can only be 1D, 2D or 3D '
                        write(*,*)'       checking D1D,D2D,D3D     '
                        stop
                endif
        else
                write(*,*)'error: please set dimension in parameters'
                stop
        endif

        ! check fluid
        if (MHD) then
                if(isothermal) then
                        write(*,*)'error: currently, we do not support'
                        write(*,*)'       isothermal MHD'
                        stop
                endif
                if(spherical) then
                        write(*,*)'error: currently, we do not support'
                        write(*,*)'       MHD in spherical coordinate'
                        stop
                endif
                if(cylindrical) then
                        write(*,*)'error: currently, we do not support'
                        write(*,*)'       MHD in cylindrical coordinate'
                        stop
                endif
        endif

        ! check geometry
        if(cartesian)then
        elseif(spherical)then
                if (.not. D1D) then
                        write(*,*)'error: currrently, spherical coord'
                        write(*,*)'       only support 1D hydro'
                        stop
                endif
                if (MHD) then
                        write(*,*)'error: currrently, spherical coord'
                        write(*,*)'       only support 1D hydro'
                        stop
                endif
        elseif(cylindrical)then
                write(*,*)'error: this version of Antares does not'
                write(*,*)'       support cylindrical coord'
                stop
        else
                write(*,*)'error: please set a geometry'
                stop
        endif
        end
