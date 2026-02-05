module grid_module
    implicit none
    private
    public :: grid_type, initialize_grid, redistribute_cells, check_critical, &
              get_sum, write_grid_state, add_grain, write_initial_state, write_grid_diff

    type :: grid_type
        integer, allocatable :: cells(:,:,:)
        integer, allocatable :: prev_cells(:,:,:)
        integer :: nx, ny, nz
        integer :: iteration
    end type grid_type

contains

    subroutine initialize_grid(grid, nx, ny, nz, init_type)
        type(grid_type), intent(out) :: grid
        integer, intent(in) :: nx, ny, nz
        character(len=*), intent(in) :: init_type
        integer :: i, j, k
        real :: rand_val

        grid%nx = nx
        grid%ny = ny
        grid%nz = nz
        grid%iteration = 0

        allocate(grid%cells(nx, ny, nz))
        allocate(grid%prev_cells(nx, ny, nz))

        select case(trim(init_type))
        case('blank')
            grid%cells = 0

        case('random')
            call random_seed()
            do i = 1, nx
                do j = 1, ny
                    do k = 1, nz
                        call random_number(rand_val)
                        grid%cells(i,j,k) = int(rand_val * 5)
                    end do
                end do
            end do

        case default
            print *, 'Unknown initialization type'
            stop
        end select

        grid%prev_cells = grid%cells

    end subroutine initialize_grid


    subroutine redistribute_cells(grid, unit)
        type(grid_type), intent(inout) :: grid
        integer, intent(in) :: unit
        integer :: i, j, k
        logical :: redistributed

        grid%iteration = grid%iteration + 1

        do while (.true.)
            redistributed = .false.
            do i = 1, grid%nx
                do j = 1, grid%ny
                    do k = 1, grid%nz
                        if (grid%cells(i,j,k) >= 6) then
                            grid%cells(i,j,k) = grid%cells(i,j,k) - 6
                            if (i > 1) grid%cells(i-1,j,k) = grid%cells(i-1,j,k) + 1
                            if (i < grid%nx) grid%cells(i+1,j,k) = grid%cells(i+1,j,k) + 1
                            if (j > 1) grid%cells(i,j-1,k) = grid%cells(i,j-1,k) + 1
                            if (j < grid%ny) grid%cells(i,j+1,k) = grid%cells(i,j+1,k) + 1
                            if (k > 1) grid%cells(i,j,k-1) = grid%cells(i,j,k-1) + 1
                            if (k < grid%nz) grid%cells(i,j,k+1) = grid%cells(i,j,k+1) + 1
                            redistributed = .true.
                        end if
                    end do
                end do
            end do
            call write_grid_diff(grid, unit)
            if (.not. redistributed) exit
        end do

    end subroutine redistribute_cells


    function check_critical(grid) result(has_critical)
        type(grid_type), intent(in) :: grid
        logical :: has_critical
        integer :: i, j, k

        has_critical = .false.
        do i = 1, grid%nx
            do j = 1, grid%ny
                do k = 1, grid%nz
                    if (grid%cells(i,j,k) >= 6) then
                        has_critical = .true.
                        return
                    end if
                end do
            end do
        end do
    end function check_critical


    subroutine add_grain(grid, x, y, z)
        type(grid_type), intent(inout) :: grid
        integer, intent(in) :: x, y, z

        ! TODO: Add 1 grain to the cell at position (x, y, z)
        ! Include bounds checking to ensure x, y, z are within grid

        if (x < 1 .or. x > grid%nx .or. &
            y < 1 .or. y > grid%ny .or. &
            z < 1 .or. z > grid%nz) then
            print *, 'Error: Coordinates out of bounds'
            return
        end if
        grid%cells(x, y, z) = grid%cells(x, y, z) + 1

    end subroutine add_grain


    function get_sum(grid) result(total)
        type(grid_type), intent(in) :: grid
        integer :: total
        integer :: i, j, k

        total = 0
        do i = 1, grid%nx
            do j = 1, grid%ny
                do k = 1, grid%nz
                    total = total + grid%cells(i,j,k)
                end do
            end do
        end do

    end function get_sum


    subroutine write_grid_state(grid, unit)
        type(grid_type), intent(in) :: grid
        integer, intent(in) :: unit
        integer :: i, j, k

        write(unit,'(A,I0)') '#', grid%iteration
        do i = 1, grid%nx
            do j = 1, grid%ny
                do k = 1, grid%nz
                    if (grid%cells(i,j,k) > 0) then
                        write(unit,'(I0,A,I0,A,I0,A,I0)') i, ',', j, ',', k, ',', grid%cells(i,j,k)
                    end if
                end do
            end do
        end do

    end subroutine write_grid_state


    subroutine write_initial_state(grid, unit)
        type(grid_type), intent(in) :: grid
        integer, intent(in) :: unit
        integer :: i, j, k

        write(unit,'(A)') '#INIT'
        do i = 1, grid%nx
            do j = 1, grid%ny
                do k = 1, grid%nz
                    if (grid%cells(i,j,k) > 0) then
                        write(unit,'(I0,A,I0,A,I0,A,I0)') i, ',', j, ',', k, ',', grid%cells(i,j,k)
                    end if
                end do
            end do
        end do

    end subroutine write_initial_state


    subroutine write_grid_diff(grid, unit)
        type(grid_type), intent(inout) :: grid
        integer, intent(in) :: unit
        integer :: i, j, k
        logical :: has_diff

        has_diff = .false.

        ! Check if there are any differences
        do i = 1, grid%nx
            do j = 1, grid%ny
                do k = 1, grid%nz
                    if (grid%cells(i,j,k) /= grid%prev_cells(i,j,k)) then
                        has_diff = .true.
                        exit
                    end if
                end do
                if (has_diff) exit
            end do
            if (has_diff) exit
        end do

        ! Write header only if there are differences
        if (has_diff) then
            write(unit,'(A,I0)') '#D', grid%iteration
            do i = 1, grid%nx
                do j = 1, grid%ny
                    do k = 1, grid%nz
                        if (grid%cells(i,j,k) /= grid%prev_cells(i,j,k)) then
                            write(unit,'(I0,A,I0,A,I0,A,I0)') i, ',', j, ',', k, ',', grid%cells(i,j,k)
                        end if
                    end do
                end do
            end do
        end if

        ! Update previous state
        grid%prev_cells = grid%cells

    end subroutine write_grid_diff

end module grid_module
