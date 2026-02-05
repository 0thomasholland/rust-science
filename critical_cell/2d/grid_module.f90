module grid_module
    implicit none
    private
    public :: grid_type, initialize_grid, redistribute_cells, check_critical, &
              get_sum, write_grid_state, add_grain

    type :: grid_type
        integer, allocatable :: cells(:,:)
        integer :: nx, ny
        integer :: iteration
    end type grid_type

contains

    subroutine initialize_grid(grid, nx, ny, init_type)
        type(grid_type), intent(out) :: grid
        integer, intent(in) :: nx, ny
        character(len=*), intent(in) :: init_type
        integer :: i, j
        real :: rand_val

        grid%nx = nx
        grid%ny = ny
        grid%iteration = 0

        allocate(grid%cells(nx, ny))

        select case(trim(init_type))
        case('blank')
            grid%cells = 0

        case('random')
            call random_seed()
            do i = 1, nx
                do j = 1, ny
                    call random_number(rand_val)
                    grid%cells(i,j) = int(rand_val * 3)
                end do
            end do

        case default
            print *, 'Unknown initialization type'
            stop
        end select

    end subroutine initialize_grid


    subroutine redistribute_cells(grid)
        type(grid_type), intent(inout) :: grid
        integer :: i, j
        logical :: redistributed

        grid%iteration = grid%iteration + 1

        do while (.true.)
            redistributed = .false.
            do i = 1, grid%nx
                do j = 1, grid%ny
                    if (grid%cells(i,j) >= 4) then
                        grid%cells(i,j) = grid%cells(i,j) - 4
                        if (i > 1) grid%cells(i-1,j) = grid%cells(i-1,j) + 1
                        if (i < grid%nx) grid%cells(i+1,j) = grid%cells(i+1,j) + 1
                        if (j > 1) grid%cells(i,j-1) = grid%cells(i,j-1) + 1
                        if (j < grid%ny) grid%cells(i,j+1) = grid%cells(i,j+1) + 1
                        redistributed = .true.
                    end if
                end do
            end do
            if (.not. redistributed) exit
        end do

    end subroutine redistribute_cells


    function check_critical(grid) result(has_critical)
        type(grid_type), intent(in) :: grid
        logical :: has_critical
        integer :: i, j

        has_critical = .false.
        do i = 1, grid%nx
            do j = 1, grid%ny
                if (grid%cells(i,j) >= 4) then
                    has_critical = .true.
                    return
                end if
            end do
        end do
    end function check_critical


    subroutine add_grain(grid, x, y)
        type(grid_type), intent(inout) :: grid
        integer, intent(in) :: x, y
        if (x < 1 .or. x > grid%nx .or. &
            y < 1 .or. y > grid%ny) then
            print *, 'Error: Coordinates out of bounds'
            return
        end if
        grid%cells(x, y) = grid%cells(x, y) + 1

    end subroutine add_grain


    function get_sum(grid) result(total)
        type(grid_type), intent(in) :: grid
        integer :: total
        integer :: i, j

        total = 0
        do i = 1, grid%nx
            do j = 1, grid%ny
                total = total + grid%cells(i,j)
            end do
        end do

    end function get_sum


    subroutine write_grid_state(grid, unit)
        type(grid_type), intent(in) :: grid
        integer, intent(in) :: unit
        integer :: i, j

        write(unit,*) 'Iteration:', grid%iteration
        do i = 1, grid%nx
            do j = 1, grid%ny
                write(unit,*) i, j, grid%cells(i,j)
            end do
        end do

    end subroutine write_grid_state

end module grid_module
