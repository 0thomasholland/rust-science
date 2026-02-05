program critical_cellular_automaton
    use grid_module
    implicit none

    type(grid_type) :: grid
    integer :: iteration, max_iterations
    integer :: random_x, random_y, random_z
    integer :: output_unit
    character(len=20) :: init_type
    character(len=256) :: output_file
    real :: rand_x, rand_y, rand_z

    integer :: nx, ny, nz

    nx = 20
    ny = 20
    nz = 20
    max_iterations = 1000
    init_type = 'random'
    output_file = 'output.csv'

    call initialize_grid(grid, nx, ny, nz, init_type)

    open(newunit=output_unit, file=trim(output_file), status='replace', action='write')

    print *, 'Starting simulation...'
    print *, 'Grid size: ', nx, ' x ', ny, ' x ', nz
    print *, 'Initialization type: ', trim(init_type)
    print *, 'Max iterations: ', max_iterations

    ! Write initial state
    call write_initial_state(grid, output_unit)

    do iteration = 1, max_iterations

        call random_number(rand_x)
        call random_number(rand_y)
        call random_number(rand_z)
        random_x = int(rand_x * nx) + 1
        random_y = int(rand_y * ny) + 1
        random_z = int(rand_z * nz) + 1
        call add_grain(grid, random_x, random_y, random_z)

        grid%iteration = iteration
        call write_grid_diff(grid, output_unit)

        do while (check_critical(grid))
            call redistribute_cells(grid, output_unit)
        end do

        if (mod(iteration, 10) == 0) then
            print *, 'Iteration: ', iteration, ' Total grains: ', get_sum(grid)
        end if

    end do

    close(output_unit)

    print *, 'Simulation complete!'
    print *, 'Final total grains: ', get_sum(grid)
    print *, 'Output written to: ', trim(output_file)

contains

end program critical_cellular_automaton
