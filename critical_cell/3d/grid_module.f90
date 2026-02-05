module grid_module
    implicit none
    private
    public :: grid_type, initialize_grid, redistribute_cells, check_critical, &
              get_sum, write_grid_state, add_grain

    type :: grid_type
        integer, allocatable :: cells(:,:,:)
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

    end subroutine initialize_grid


    subroutine redistribute_cells(grid)
        type(grid_type), intent(inout) :: grid
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

        write(unit,*) 'Iteration:', grid%iteration
        do i = 1, grid%nx
            do j = 1, grid%ny
                do k = 1, grid%nz
                    write(unit,*) i, j, k, grid%cells(i,j,k)
                end do
            end do
        end do

    end subroutine write_grid_state

end module grid_module
