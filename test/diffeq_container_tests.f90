program diffeq_container_tests
    use iso_fortran_env, only : real64
    use diffeq_containers
    implicit none

    character(len=40) :: test_name
    type(ode_results) :: results
    real(real64) :: point
    real(real64), allocatable, dimension(:) :: values, grid, vector
    real(real64), allocatable, dimension(:,:) :: matrix

    call get_command_argument(1, test_name)
    if (len_trim(test_name) == 0) then
        call test_construction()
        call test_state_access()
        call test_capacity()
        call test_push()
        call test_to_matrix()
        call test_clear_all()
    else
        select case (trim(test_name))
        case ("init_zero_states")
            results = ode_results(0)
        case ("init_negative_states")
            results = ode_results(-1)
        case ("init_zero_capacity")
            results = ode_results(2, 0)
        case ("init_negative_capacity")
            results = ode_results(2, -1)
        case ("capacity_uninitialized")
            call results%set_capacity(1)
        case ("point_uninitialized")
            point = results%solution_point(1)
        case ("selected_uninitialized")
            values = results%selected_state(1)
        case ("whole_uninitialized")
            matrix = results%state()
        case ("vector_uninitialized")
            vector = results%state(1)
        case ("points_uninitialized")
            grid = results%solution_points()
        case ("matrix_uninitialized")
            matrix = results%to_matrix()
        case ("point_empty", "selected_empty")
            results = ode_results(2)
            if (trim(test_name) == "point_empty") then
                point = results%solution_point(1)
            else
                values = results%selected_state(1)
            end if
        case default
            results = ode_results(2, 2)
            call results%push(1.0d0, [2.0d0, 3.0d0])
            select case (trim(test_name))
            case ("capacity_zero")
                call results%set_capacity(0)
            case ("capacity_negative")
                call results%set_capacity(-1)
            case ("push_size")
                call results%push(2.0d0, [4.0d0])
            case ("point_low")
                point = results%solution_point(0)
            case ("point_high")
                point = results%solution_point(2)
            case ("selected_low")
                values = results%selected_state(0)
            case ("selected_high")
                values = results%selected_state(2)
            case ("vector_low")
                vector = results%state(0)
            case ("vector_high")
                vector = results%state(3)
            case ("point_cleared")
                call results%clear_all()
                point = results%solution_point(1)
            case ("selected_cleared")
                call results%clear_all()
                values = results%selected_state(1)
            case ("vector_cleared")
                call results%clear_all()
                vector = results%state(1)
            case default
                error stop "Unknown container test"
            end select
        end select
        stop 0
    end if

contains
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not.condition) error stop message
    end subroutine

    subroutine test_construction()
        type(ode_results) :: container
        real(real64), allocatable, dimension(:) :: empty_grid, empty_vector
        real(real64), allocatable, dimension(:,:) :: empty_states

        call check(container%count() == 0, "Default count")
        call check(container%get_capacity() == 0, "Default capacity")
        call check(container%state_variable_count() == 0, "Default state count")

        container = ode_results(2)
        call check(container%count() == 0, "Constructed count")
        call check(container%get_capacity() == 1000, "Default buffer capacity")
        call check(container%state_variable_count() == 2, "Constructed state count")
        empty_states = container%state()
        empty_vector = container%state(2)
        empty_grid = container%solution_points()
        call check(allocated(empty_states), "Empty state allocated")
        call check(all(shape(empty_states) == [0, 2]), "Empty state shape")
        call check(allocated(empty_vector), "Empty state vector allocated")
        call check(size(empty_vector) == 0, "Empty state vector size")
        call check(allocated(empty_grid), "Empty solution points allocated")
        call check(size(empty_grid) == 0, "Empty solution points size")

        container = ode_results(1, 1)
        call check(container%count() == 0, "Explicit empty count")
        call check(container%get_capacity() == 1, "Minimum capacity")
        call check(container%state_variable_count() == 1, "Single state count")
    end subroutine

    subroutine test_state_access()
        type(ode_results) :: container
        real(real64) :: xvalue
        real(real64), allocatable, dimension(:) :: yvalue, xvalues, yvector
        real(real64), allocatable, dimension(:,:) :: yvalues
        real(real64), dimension(3,2) :: expected
        integer :: index

        container = ode_results(2, 5)
        expected(:,1) = [10.0d0, 20.0d0, 30.0d0]
        expected(:,2) = [-1.0d0, -2.0d0, -3.0d0]
        do index = 1, 3
            call container%push(real(index, real64), expected(index,:))
        end do
        call check(container%count() == 3, "Push count")
        call check(container%get_capacity() == 5, "Push within capacity")

        do index = 1, 3
            xvalue = container%solution_point(index)
            yvalue = container%selected_state(index)
            call check(size(yvalue) == 2, "Selected state size")
            call check(xvalue == real(index, real64), "Indexed independent variable")
            call check(all(yvalue == expected(index,:)), "Indexed state")
        end do
        xvalues = container%solution_points()
        yvalues = container%state()
        call check(size(xvalues) == 3, "Solution points exclude spare capacity")
        call check(all(shape(yvalues) == [3, 2]), "Whole state excludes spare capacity")
        call check(all(xvalues == [1.0d0, 2.0d0, 3.0d0]), "Whole independent variable")
        call check(all(yvalues == expected), "Whole state")
        do index = 1, 2
            yvector = container%state(index)
            call check(size(yvector) == 3, "State vector excludes spare capacity")
            call check(all(yvector == expected(:,index)), "State vector with spare capacity")
        end do

        yvalue = 0.0d0
        yvalue = container%selected_state(3)
        call check(all(yvalue == expected(3,:)), "Selected state is a copy")
        xvalues = 0.0d0
        xvalues = container%solution_points()
        call check(all(xvalues == [1.0d0, 2.0d0, 3.0d0]), "Solution points are a copy")
        yvalues = 0.0d0
        yvalues = container%state()
        call check(all(yvalues == expected), "Whole state is a copy")
        yvector = 0.0d0
        yvector = container%state(2)
        call check(all(yvector == expected(:,2)), "State vector is a copy")
        call container%clear_all()
        call check(all(yvalues == expected), "Clear preserves returned state copy")
        call check(all(xvalues == [1.0d0, 2.0d0, 3.0d0]), "Clear preserves returned points copy")
    end subroutine

    subroutine test_capacity()
        type(ode_results) :: container
        real(real64) :: xvalue
        real(real64), dimension(2) :: yvalue
        integer :: index

        container = ode_results(2, 3)
        call container%set_capacity(1)
        call check(container%count() == 0, "Empty resize count")
        call check(container%get_capacity() == 1, "Empty resize capacity")
        call container%set_capacity(3)
        do index = 1, 3
            call container%push(real(index, real64), &
                [real(index * 10, real64), real(-index, real64)])
        end do
        call container%set_capacity(6)
        call check(container%get_capacity() == 6, "Grow capacity")
        call check(container%count() == 3, "Grow count")
        call check(container%state_variable_count() == 2, "Resize state count")
        call container%set_capacity(6)
        call check(container%get_capacity() == 6, "Unchanged capacity")
        call check(container%count() == 3, "Unchanged count")
        call container%set_capacity(4)
        call check(container%get_capacity() == 4, "Shrink spare capacity")
        call check(container%count() == 3, "Shrink without truncation")
        do index = 1, 3
            xvalue = container%solution_point(index)
            yvalue = container%selected_state(index)
            call check(xvalue == real(index, real64), "Resize preserves independent variable")
            call check(all(yvalue == [real(index * 10, real64), real(-index, real64)]), &
                "Resize preserves state")
        end do
        call container%set_capacity(3)
        call check(container%count() == 3, "Shrink to count")
        call check(container%get_capacity() == 3, "Exact count capacity")
        call container%set_capacity(1)
        call check(container%get_capacity() == 1, "Truncated capacity")
        call check(container%count() == 1, "Truncated count")
        xvalue = container%solution_point(1)
        yvalue = container%selected_state(1)
        call check(xvalue == 1.0d0 .and. all(yvalue == [10.0d0, -1.0d0]), &
            "Truncation keeps first point")
        call container%push(4.0d0, [40.0d0, -4.0d0])
        call check(container%count() == 2, "Push after truncation")
        xvalue = container%solution_point(2)
        yvalue = container%selected_state(2)
        call check(xvalue == 4.0d0 .and. all(yvalue == [40.0d0, -4.0d0]), &
            "Push replaces discarded point")
    end subroutine

    subroutine test_clear_all()
        type(ode_results) :: container
        real(real64), allocatable, dimension(:,:) :: output
        real(real64), allocatable, dimension(:) :: empty_grid
        real(real64), allocatable, dimension(:,:) :: empty_states

        call container%clear_all()
        call check(container%count() == 0, "Clear default count")
        call check(container%get_capacity() == 0, "Clear default capacity")
        call check(container%state_variable_count() == 0, "Clear default state count")

        container = ode_results(2, 3)
        call container%clear_all()
        call check(container%count() == 0, "Clear empty count")
        call check(container%get_capacity() == 3, "Clear empty preserves capacity")
        call check(container%state_variable_count() == 0, "Clear empty state count")
        call container%push(1.0d0, [2.0d0, 3.0d0])
        call container%push(4.0d0, [5.0d0, 6.0d0])
        call container%clear_all()
        call container%clear_all()
        call check(container%count() == 0, "Repeated clear count")
        call check(container%get_capacity() == 3, "Clear populated preserves capacity")
        call check(container%state_variable_count() == 0, "Clear populated state count")
        output = container%to_matrix()
        call check(all(shape(output) == [0, 1]), "Cleared matrix shape")
        empty_grid = container%solution_points()
        empty_states = container%state()
        call check(size(empty_grid) == 0, "Cleared solution points size")
        call check(all(shape(empty_states) == [0, 0]), "Cleared state shape")

        call container%push(7.0d0, [8.0d0, 9.0d0])
        output = container%to_matrix()
        call check(all(shape(output) == [1, 3]), "Reuse cleared matrix shape")
        call check(all(output(1,:) == [7.0d0, 8.0d0, 9.0d0]), "Clear removes old points")
        call check(container%get_capacity() == 3, "Reuse retains capacity")

        call container%clear_all()
        call container%set_capacity(4)
        call container%push(10.0d0, [11.0d0])
        output = container%to_matrix()
        call check(container%get_capacity() == 4, "Resize cleared capacity")
        call check(container%state_variable_count() == 1, "Reuse with fewer variables")
        call check(all(shape(output) == [1, 2]), "Fewer variables matrix shape")
        call check(all(output(1,:) == [10.0d0, 11.0d0]), "Fewer variables matrix values")

        call container%clear_all()
        call container%push(12.0d0, [13.0d0, 14.0d0, 15.0d0])
        output = container%to_matrix()
        call check(container%get_capacity() == 4, "New vector size preserves capacity")
        call check(container%state_variable_count() == 3, "Reuse with more variables")
        call check(all(shape(output) == [1, 4]), "More variables matrix shape")
        call check(all(output(1,:) == [12.0d0, 13.0d0, 14.0d0, 15.0d0]), "More variables values")
    end subroutine

    subroutine test_to_matrix()
        type(ode_results) :: container
        real(real64), allocatable, dimension(:,:) :: output, snapshot
        real(real64), dimension(2,3) :: expected

        container = ode_results(2, 5)
        output = container%to_matrix()
        call check(allocated(output), "Empty matrix allocated")
        call check(all(shape(output) == [0, 3]), "Empty matrix shape")

        call container%push(1.0d0, [10.0d0, -1.0d0])
        call container%push(2.0d0, [20.0d0, -2.0d0])
        expected(:,1) = [1.0d0, 2.0d0]
        expected(:,2) = [10.0d0, 20.0d0]
        expected(:,3) = [-1.0d0, -2.0d0]
        output = container%to_matrix()
        call check(all(shape(output) == [2, 3]), "Matrix excludes spare capacity")
        call check(all(output == expected), "Matrix column order and values")
        call check(container%count() == 2, "Conversion preserves count")
        call check(container%get_capacity() == 5, "Conversion preserves capacity")
        call check(container%state_variable_count() == 2, "Conversion preserves state count")

        output = 0.0d0
        snapshot = container%to_matrix()
        call check(all(snapshot == expected), "Matrix is an independent copy")
        call container%set_capacity(1)
        output = container%to_matrix()
        call check(all(shape(output) == [1, 3]), "Truncated matrix shape")
        call check(all(output(1,:) == expected(1,:)), "Truncated matrix values")
        call container%push(3.0d0, [30.0d0, -3.0d0])
        output = container%to_matrix()
        call check(all(shape(output) == [2, 3]), "Matrix shape after growth")
        call check(all(output(1,:) == expected(1,:)), "Matrix retains first row after growth")
        call check(all(output(2,:) == [3.0d0, 30.0d0, -3.0d0]), "Matrix new row after growth")
        call check(all(snapshot == expected), "Container changes preserve matrix copy")
    end subroutine

    subroutine test_push()
        type(ode_results) :: container, small
        real(real64) :: xvalue
        real(real64), dimension(2) :: yvalue
        integer :: index

        do index = 1, 2001
            call container%push(real(index, real64), &
                [real(index * 2, real64), real(-index, real64)])
        end do
        call check(container%state_variable_count() == 2, "Lazy initialization")
        call check(container%count() == 2001, "Repeated growth count")
        call check(container%get_capacity() >= container%count(), "Repeated growth capacity")
        do index = 1, container%count()
            xvalue = container%solution_point(index)
            yvalue = container%selected_state(index)
            call check(xvalue == real(index, real64), "Growth preserves independent variable")
            call check(all(yvalue == [real(index * 2, real64), real(-index, real64)]), &
                "Growth preserves state")
        end do

        small = ode_results(2, 1)
        call small%push(1.0d0, [2.0d0, 3.0d0])
        call small%push(4.0d0, [5.0d0, 6.0d0])
        call check(small%count() == 2, "Explicit capacity growth count")
        call check(small%get_capacity() >= 2, "Explicit capacity growth")
        xvalue = small%solution_point(1)
        yvalue = small%selected_state(1)
        call check(xvalue == 1.0d0 .and. all(yvalue == [2.0d0, 3.0d0]), "Growth first point")
        xvalue = small%solution_point(2)
        yvalue = small%selected_state(2)
        call check(xvalue == 4.0d0 .and. all(yvalue == [5.0d0, 6.0d0]), "Growth last point")
    end subroutine
end program