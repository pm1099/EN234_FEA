subroutine user_print(n_steps)
  use Types
  use ParamIO
  use Globals, only : TIME, DTIME
  use Mesh
  use Printparameters, only : n_user_print_files                  ! No. files specified by the user
  use Printparameters, only : n_user_print_parameters             ! No. user supplied print parameters
  use Printparameters, only : user_print_units                    ! Unit numbers
  use Printparameters, only : user_print_parameters               ! List of user supplied parameters
  use User_Subroutine_Storage, only : length_state_variable_array ! Max no. state variables on any element
  implicit none
  
  integer, intent(in) :: n_steps                                 ! Current step number
!  real (prec) :: j_integral_value
  integer ::  lmn
  integer ::  status
  integer ::  n_state_vars_per_intpt                                         ! No. state variables per integration point
  real (prec) ::   vol_averaged_strain(6)                                    ! Volume averaged strain in an element
  real (prec) ::   vol_averaged_stress(6)
  real (prec), allocatable ::   vol_averaged_state_variables(:)              ! Volume averaged state variables in an element



!
!  Use this file to process or print time histories of the solution, or to print a non-standard mesh.
!
!  As an example, this subroutine computes the volume averaged infinitesimal strain and the volume average
!  element state variables (if they exist) in an element.   The element is specified by the user.
!  The first six state variables (which are usually the stresses) are printed along with the strains.
!
!

   allocate(vol_averaged_state_variables(length_state_variable_array), stat=status)

   if (status/=0) then
      write(IOW,*) ' Error in subroutine user_print'
     write(IOW,*) ' Unable to allocate memory for state variables '
      stop
   endif

   lmn = int(user_print_parameters(1))     ! The element number
   !call compute_J_integral (J_integral_value)
   !write (6,*) J_integral_value

   call compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_stress,vol_averaged_stress,&
   vol_averaged_state_variables,length_state_variable_array,n_state_vars_per_intpt)

   write(user_print_units(1),'(2(D15.5))') vol_averaged_stress(1),vol_averaged_strain(1)

!    if (TIME<1.d-12) then
!      if (n_state_vars_per_intpt<6) then
!        write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23'
!      else
!         write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23'
!      endif
!    endif
!
!   if (n_state_vars_per_intpt<6) then
!      write(user_print_units(1),'(7(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6)
!   else
!      vol_averaged_state_variables(1:3) = vol_averaged_state_variables(1:3) + vol_averaged_state_variables(7)
!     write(user_print_units(1),'(13(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6),vol_averaged_stress(1:6)
!   endif



end subroutine user_print
!
subroutine compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_stress,vol_averaged_state_vars,&
length_output_array,n_state_vars_per_intpt)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use User_Subroutine_Storage
    use Element_Utilities, only : Ns => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent ( in )      :: lmn                                          ! Element number
    integer, intent ( in )      :: length_output_array

    real (prec), intent( out )  ::  vol_averaged_strain(6)
    real (prec), intent( out )  ::  vol_averaged_stress(6)
    real (prec), intent( out )  ::  vol_averaged_state_vars(length_output_array)

    integer, intent( out )      :: n_state_vars_per_intpt

    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element


    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    integer      :: n_points,kint,i
    integer      :: n_coords, n_dof
    integer      :: iof
    integer      :: status

    real (prec)  ::  el_vol
    real (prec), allocatable  ::  B(:,:)               ! strain = B*(dof_total+dof_increment)
    real (prec), allocatable  ::  Bcorr(:,:)
    real (prec), allocatable  ::  Bbar(:,:)
    real (prec)  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  dstrain(6)                        ! Strain increment vector
    real (prec)  ::  strain_total(6),strain_dev(6)     ! strain_total=strain+dstrain
    real (prec)  ::  stress(6)
    real (prec)  ::  sigma_0,epsilon_0,n,K             ! Element Properties
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  epsilon_e,sigma_e                 ! Equivalent Strain and equivalent stress
    real (prec)  :: D(6,6)
    !
    !  Allocate memory to store element data.
    !  The variables specifying the size of the arrays are stored in the module user_subroutine_storage
    !  They are initialized when the input file is read, and specify the size of the arrays required to store data
    !  for any element in the mesh.  Some elements may require less storage.

    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(6,length_dof_array), stat=status)
    allocate(Bcorr(6,length_dof_array), stat=status)
    allocate(Bbar(6,length_dof_array), stat=status)

    if (status/=0) then
       write(IOW,*) ' Error in subroutine compute_volume_average_3D'
       write(IOW,*) ' Unable to allocate memory for element variables '
       stop
    endif
    !
    ! Extract element and node data from global storage (see module Mesh.f90 for the source code for these subroutines)

    call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

    do i = 1, n_nodes
        iof = 3*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:3,i),n_dof, &
                                                 dof_increment(iof:iof+2),dof_total(iof:iof+2))
    end do

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    vol_averaged_strain = 0.d0
    vol_averaged_stress = 0.d0
!    vol_averaged_state_vars = 0.d0
    el_vol = 0.d0


    n_state_vars_per_intpt = 0

    if (n_state_vars_per_intpt>size(vol_averaged_state_vars)) then
       write(IOW,*) ' Error detected in subroutine compute_element_volume_average_3d '
       write(IOW,*) ' The element contains ',n_state_vars_per_intpt
       write(IOW,*) ' but the array storing averaged state variables has length ',size(vol_averaged_state_vars)
       stop
    endif

    !Element Properties
    sigma_0 = element_properties(1)
    epsilon_0 = element_properties(2)
    n=element_properties(3)
    K=element_properties(4)

    ! Finding element volume and B-bar matrix. This needs to be done outside the main integration points loop as the main loop requires the value of B-bar matrix at each integration point.
    el_vol=0.d0
    dNbardx=0.d0

    do kint=1,n_points
        ! Find dNdx at each integration point
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,Ns,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dNbardx(1:n_nodes,1:3)=dNbardx(1:n_nodes,1:3)+dNdx(1:n_nodes,1:3)*w(kint)*determinant
        el_vol=el_vol+w(kint)*determinant
    end do
    dNbardx=dNbardx/el_vol
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,Ns,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)

        iof = n_state_vars_per_intpt*(kint-1)+1
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        ! Finding Bcorr matrix
        Bcorr=0.d0
        Bcorr(1,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1)
        Bcorr(1,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2)
        Bcorr(1,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3)
        Bcorr(2,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1)
        Bcorr(2,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2)
        Bcorr(2,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3)
        Bcorr(3,1:3*n_nodes-2:3)=dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1)
        Bcorr(3,2:3*n_nodes-1:3)=dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2)
        Bcorr(3,3:3*n_nodes:3)=dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3)

        ! Find Bbar matrix
        Bbar=0.d0
        Bbar=B+(1/3.d0)*Bcorr

        ! Finding strain.
        strain = matmul(Bbar(1:6,1:3*n_nodes),dof_total(1:3*n_nodes))
        dstrain = matmul(Bbar(1:6,1:3*n_nodes),dof_increment(1:3*n_nodes))
        strain_total=strain+dstrain

        call tanstiff(D,stress,strain_total,sigma_0,epsilon_0,n,K)

        vol_averaged_strain(1:6) = vol_averaged_strain(1:6) + (strain_total(1:6))*w(kint)*determinant
        vol_averaged_stress(1:6) = vol_averaged_stress(1:6) + stress(1:6)*w(kint)*determinant

!        if (n_state_vars_per_intpt>0) then
!           vol_averaged_state_vars(1:n_state_vars_per_intpt) = vol_averaged_state_vars(1:n_state_vars_per_intpt) &
!                              + updated_state_variables(iof:iof+n_state_vars_per_intpt-1)*w(kint)*determinant
!        endif

    end do

    vol_averaged_strain = vol_averaged_strain/el_vol
    vol_averaged_stress = vol_averaged_stress/el_vol
 !   vol_averaged_state_vars = vol_averaged_state_vars/el_vol

    deallocate(node_list)
    deallocate(element_properties)
   deallocate(initial_state_variables)
   deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)
    deallocate(Bcorr)
    deallocate(Bbar)

    return




end subroutine compute_element_volume_average_3D

subroutine compute_J_integral(J_integral_value)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use Mesh, only : zone,zone_list
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    real (prec), intent( out )  ::  J_integral_value


    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element
    integer    :: n_coords                                     ! No. coords for a node
    integer    :: n_dof                                        ! No. DOFs for a node
    integer      :: status
    integer      :: iof
    integer      :: lmn               ! Element number
    integer      :: lmn_start,lmn_end ! First and last crack tip element
    integer      :: i                 ! Loop counter
    integer        ::n_points, kint           !Number of integration points and loop counter for integration points

    real (prec)     ::xint(2)                 ! Coordinates of integration points
    real (prec)     ::qderiv(2)               ! derivatives of q
    real (prec)     ::u12                     ! u12=du1/du2
    real (prec)     ::r                       ! r=sqrt(xint(1)^2+xint(2)^2)
    real (prec)     ::sed                       ! Strain energy density

    real (prec)  ::  strain(3)            ! Strain vector contains [e11, e22, 2e12]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  E, xnu, d33, d11, d12             ! Material properties

!   The arrays below have to be given dimensions large enough to store the data. It doesnt matter if they are too large.

    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    real (prec), allocatable  ::  B(:,:)                                   ! strain = B*(dof_total+dof_increment)
    !
    !
    !  The variables specifying the sizes of the arrays listed below are determined while reading the input file
    !  They specify the array dimensions required to store the relevant variables for _any_ element or node in the mesh
    !  The actual data will vary depending on the element or node selected
    !
    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(2,length_coord_array/2), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(3,length_dof_array), stat=status)

  !  Write your code to calculate the J integral here
    J_integral_value=0.d0
  !  You will need to loop over the crack tip elements, and sum the contribution to the J integral from each element.
  !
  !  You can access the first and last crack tip element using
     lmn_start = zone_list(2)%start_element
     lmn_end = zone_list(2)%end_element

    ! Loop over all the crack tip elements to find the J-integral value.
     do lmn=lmn_start,lmn_end

        ! Extract element data
        call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)
        if (n_nodes == 6) n_points = 3  ! 6-noded triangular elements
        if (n_nodes == 8) n_points = 9

        ! Extract node data
        do i = 1, n_nodes
        iof = 2*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:2,i),n_dof, &
                                                 dof_increment(iof:iof+1),dof_total(iof:iof+1))
        end do

        call initialize_integration_points(n_points, n_nodes, xi, w)

        D = 0.d0
        E = element_properties(1)
        xnu = element_properties(2)
        d33 = 0.5D0*E/(1+xnu)
        d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
        d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
        D(1:2,1:2) = d12
        D(1,1) = d11
        D(2,2) = d11
        D(3,3) = d33

        sed=0.d0
        ! Find the J-integral for element lmn by adding contributions from each integration point.
        do kint = 1, n_points

            call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx(1:2,1:2))
            B = 0.d0
            B(1,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
            B(2,2:2*n_nodes:2) = dNdx(1:n_nodes,2)
            B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
            B(3,2:2*n_nodes:2)   = dNdx(1:n_nodes,1)

            strain = matmul(B,dof_total+dof_increment)
            stress = matmul(D,strain)

            xint=matmul(x(1:2,1:n_nodes),N(1:n_nodes))
            if (r>0.0006) then
            qderiv=0.d0
            else
            r=sqrt((xint(1)**2.d0)+(xint(2)**2.d0))
            qderiv=-xint/(0.0006d0*r)
            endif
            u12=dot_product((dNdx(1:n_nodes,2)),dof_total(1:2*n_nodes-1:2)+dof_increment(1:2*n_nodes-1:2))
            sed=0.5d0*(stress(1)*strain(1)+stress(2)*strain(2)+2*stress(3)*strain(3))

            !J_integral_value=J_integral_value+(((stress(1)*u12+stress(3)*strain(2))*qderiv(1))+((stress(3)*u12+stress(2)*strain(2)-sed)*qderiv(2)))*determinant*w(kint)
            J_integral_value=J_integral_value+((stress(1)*u12+stress(3)*strain(2))*qderiv(1))*determinant*w(kint)&
            +((stress(3)*u12+stress(2)*strain(2)-sed)*qderiv(2))*determinant*w(kint)
        end do

     end do

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)


    return




end subroutine compute_J_integral

