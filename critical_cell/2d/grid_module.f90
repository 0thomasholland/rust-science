module grid_module
    implicit none
    private
    public :: grid_type, initialize_grid, redistribute_cells, check_critical, &
              get_sum, write_grid_state, add_grain, write_initial_state, write_grid_diff

    type :: grid_type
        integer, allocatable :: cells(:,:)
        integer, allocatable :: prev_cells(:,:)
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
        allocate(grid%prev_cells(nx, ny))

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

        grid%prev_cells = grid%cells

    end subroutine initialize_grid


    subroutine redistribute_cells(grid, unit)
        type(grid_type), intent(inout) :: grid
        integer, intent(in) :: unit
        integer, allocatable :: temp_cells(:,:)
        integer :: i, j
        logical :: redistributed

        grid%iteration = grid%iteration + 1
        allocate(temp_cells(grid%nx, grid%ny))

        do while (.true.)
            redistributed = .false.

            ! Phase 1: Process even cells (i+j is even) - no atomic contention
            temp_cells = grid%cells
            !$omp parallel do collapse(2) reduction(.or.:redistributed)
            do i = 1, grid%nx
                do j = 1, grid%ny
                    if (mod(i+j, 2) == 0 .and. temp_cells(i,j) >= 4) then
                        temp_cells(i,j) = temp_cells(i,j) - 4

                        if (i > 1) temp_cells(i-1,j) = temp_cells(i-1,j) + 1
                        if (i < grid%nx) temp_cells(i+1,j) = temp_cells(i+1,j) + 1
                        if (j > 1) temp_cells(i,j-1) = temp_cells(i,j-1) + 1
                        if (j < grid%ny) temp_cells(i,j+1) = temp_cells(i,j+1) + 1
                        redistributed = .true.
                    end if
                end do
            end do
            !$omp end parallel do
            grid%cells = temp_cells

            ! Phase 2: Process odd cells (i+j is odd) - no atomic contention
            temp_cells = grid%cells
            !$omp parallel do collapse(2) reduction(.or.:redistributed)
            do i = 1, grid%nx
                do j = 1, grid%ny
                    if (mod(i+j, 2) == 1 .and. temp_cells(i,j) >= 4) then
                        temp_cells(i,j) = temp_cells(i,j) - 4

                        if (i > 1) temp_cells(i-1,j) = temp_cells(i-1,j) + 1
                        if (i < grid%nx) temp_cells(i+1,j) = temp_cells(i+1,j) + 1
                        if (j > 1) temp_cells(i,j-1) = temp_cells(i,j-1) + 1
                        if (j < grid%ny) temp_cells(i,j+1) = temp_cells(i,j+1) + 1
                        redistributed = .true.
                    end if
                end do
            end do
            !$omp end parallel do
            grid%cells = temp_cells

            call write_grid_diff(grid, unit)
            if (.not. redistributed) exit
        end do

        deallocate(temp_cells)

    end subroutine redistribute_cells


    function check_critical(grid) result(has_critical)
        type(grid_type), intent(in) :: grid
        logical :: has_critical
        integer :: i, j

        has_critical = .false.
        !$omp parallel do collapse(2) reduction(.or.:has_critical)
        do i = 1, grid%nx
            do j = 1, grid%ny
                if (grid%cells(i,j) >= 4) then
                    has_critical = .true.
                end if
            end do
        end do
        !$omp end parallel do
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
        !$omp parallel do collapse(2) reduction(+:total)
        do i = 1, grid%nx
            do j = 1, grid%ny
                total = total + grid%cells(i,j)
            end do
        end do
        !$omp end parallel do

    end function get_sum


    subroutine write_grid_state(grid, unit)
        type(grid_type), intent(in) :: grid
        integer, intent(in) :: unit
        integer :: i, j

        write(unit,'(A,I0)') '#', grid%iteration
        do i = 1, grid%nx
            do j = 1, grid%ny
                if (grid%cells(i,j) > 0) then
                    write(unit,'(I0,A,I0,A,I0)') i, ',', j, ',', grid%cells(i,j)
                end if
            end do
        end do

    end subroutine write_grid_state


    subroutine write_initial_state(grid, unit)
        type(grid_type), intent(in) :: grid
        integer, intent(in) :: unit
        integer :: i, j

        write(unit,'(A)') '#INIT'
        do i = 1, grid%nx
            do j = 1, grid%ny
                if (grid%cells(i,j) > 0) then
                    write(unit,'(I0,A,I0,A,I0)') i, ',', j, ',', grid%cells(i,j)
                end if
            end do
        end do

    end subroutine write_initial_state


    subroutine write_grid_diff(grid, unit)
        type(grid_type), intent(inout) :: grid
        integer, intent(in) :: unit
        integer :: i, j
        logical :: has_diff

        has_diff = .false.

        ! Single pass: write diffs as we find them, header written on first diff
        do i = 1, grid%nx
            do j = 1, grid%ny
                if (grid%cells(i,j) /= grid%prev_cells(i,j)) then
                    if (.not. has_diff) then
                        write(unit,'(A,I0)') '#D', grid%iteration
                        has_diff = .true.
                    end if
                    write(unit,'(I0,A,I0,A,I0)') i, ',', j, ',', grid%cells(i,j)
                end if
            end do
        end do

        ! Update previous state
        grid%prev_cells = grid%cells

    end subroutine write_grid_diff

end module grid_module
