!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    use Globals, only: TIME,DTIME  !For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  Bcorr(6,length_dof_array)         ! Correction Matrix
    real (prec)  ::  Bbar(6,length_dof_array)          ! Bbar=B+(1/3)*Bcorr
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, D44, D11, D12              ! Material properties
    real (prec)  ::  el_vol              ! Element Volume
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44

    if (element_identifier==1001)then
    !     --  Loop over integration points
        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
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

            strain = matmul(B,dof_total)
            dstrain = matmul(B,dof_increment)

            stress = matmul(D,strain+dstrain)
            element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

            element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
                + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

        end do
    else if (element_identifier==1002) then
        ! Finding element volume and B-bar matrix. This needs to be done outside the main integration points loop as the main loop requires the value of B-bar matrix at each integration point.
        el_vol=0.d0
        dNbardx=0.d0
        do kint=1,n_points
            ! Find dNdx at each integration point
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
            dNbardx(1:n_nodes,1:3)=dNbardx(1:n_nodes,1:3)+dNdx(1:n_nodes,1:3)*w(kint)*determinant
            el_vol=el_vol+w(kint)*determinant
        end do
        dNbardx=dNbardx/el_vol

        ! Now write the main loop.
        do kint=1,n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
            ! Find the original B matrix.
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
            Bcorr(3,3:3*n_nodes:2)=dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3)

            ! Find Bbar matrix
            Bbar=0.d0
            Bbar=B+(1/3)*Bcorr

            strain = matmul(Bbar,dof_total)
            dstrain = matmul(Bbar,dof_increment)

            stress = matmul(D,strain+dstrain)
            element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(Bbar),stress)*w(kint)*determinant

            element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
                + matmul(transpose(Bbar(1:6,1:3*n_nodes)),matmul(D,Bbar(1:6,1:3*n_nodes)))*w(kint)*determinant

        end do


    endif
    return
end subroutine el_linelast_3dbasic


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_linelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
  
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
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

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      
        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
  
    return
end subroutine el_linelast_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only: dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  Bcorr(6,length_dof_array)             ! Correction Matrix
    real (prec)  ::  Bbar(6,length_dof_array)             ! Bbar=B+(1/3)*Bcorr
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: p, smises                          ! Pressure and Mises stress
    real (prec)  :: el_vol                             ! Element Volume
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
  
    if (element_identifier==1001) then
        !     --  Loop over integration points
        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
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

            strain = matmul(B,dof_total)
            dstrain = matmul(B,dof_increment)
            stress = matmul(D,strain+dstrain)
            p = sum(stress(1:3))/3.d0
            sdev = stress
            sdev(1:3) = sdev(1:3)-p
            smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
            ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
            do k = 1,n_field_variables
                if (strcmp(field_variable_names(k),'S11',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(1)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'S22',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(2)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'S33',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(3)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'S12',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(4)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'S13',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(5)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'S23',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(6)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + smises*N(1:n_nodes)*determinant*w(kint)
                endif
            end do

        end do

    else if (element_identifier==1002) then
        ! Finding element volume and B-bar matrix. This needs to be done outside the main integration points loop as the main loop requires the value of B-bar matrix at each integration point.
        el_vol=0.d0
        dNbardx=0.d0
        do kint=1,n_points
            ! Find dNdx at each integration point
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
            dNbardx(1:n_nodes,1:3)=dNbardx(1:n_nodes,1:3)+dNdx(1:n_nodes,1:3)*w(kint)*determinant
            el_vol=el_vol+w(kint)*determinant
        end do
        dNbardx=dNbardx/el_vol

        ! Now write the main loop.
        do kint=1,n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
            ! Find the original B matrix.
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
            Bbar=B+(1/3)*Bcorr

            strain = matmul(Bbar,dof_total)
            dstrain = matmul(Bbar,dof_increment)

            stress = matmul(D,strain+dstrain)
                        p = sum(stress(1:3))/3.d0
            sdev = stress
            sdev(1:3) = sdev(1:3)-p
            smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
            ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
            do k = 1,n_field_variables
                if (strcmp(field_variable_names(k),'S11',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(1)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'S22',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(2)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'S33',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(3)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'S12',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(4)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'S13',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(5)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'S23',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + stress(6)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) &
                    + smises*N(1:n_nodes)*determinant*w(kint)
                endif
            end do
        end do
    endif
    return
end subroutine fieldvars_linelast_3dbasic

! UMAT for 3D small strain viscoplasticity.

SUBROUTINE usermat_viscoplastic(STRESS,STATEV,DDSDDE,STRAN,DSTRAN,&
            PROPS,NPROPS,USTATV)

     use Types
     use ParamIO
     use Globals, only: TIME,DTIME
     use Mesh, only : node

     implicit none

     integer, intent (in)   ::  NPROPS        !7 for 3D viscoplastic

     real (prec), intent (in)        :: STRAN(6), DSTRAN(6)
     real (prec), intent (in)        :: PROPS(NPROPS)
     real (prec), intent (inout)     :: STATEV(1)

     real (prec), intent (out)  ::  STRESS(6), DDSDDE(6,6), USTATV(1)

     ! Local variables

     integer k,NIT,MAXIT,i,j                                ! Integers
     integer n,m                                            ! Integer material properties

     real(prec):: E,xnu,Y,xe,edot                   ! Material Properties
     real(prec):: EP                                ! State variables: Plastic strain
     real(prec):: AONE(6,6), ATWO(6,6)              ! AONE=dik*djl+djk*il, ATWO=dij*dkl
     real(prec):: DDEVST(6), DEVS(6)        ! Deviatoric strain increment and stress
     real(prec):: DEVSS(6), SES                 ! Elastic Predictors
     real(prec):: ERROR,TOL                           ! Tolerance on Newton-Raphson loop
     real(prec):: DEP,T,F,DFDE                          ! N-R loop terms
     real(prec):: BETA,GAMMAA, FACTOR                    ! GAMMAA
     real(prec):: ATHREE(6,6)                       ! =Sij*Skl


    ! Material Properties
    E = props(1)
    xnu = props(2)
    Y = props(3)
    xe = props(4)
    n = props(5)
    edot = props(6)
    m = props(7)

    ! State variables
    EP = STATEV(1)           ! Plane strain at nth time


    ! Define AONE and ATWO matrices
    AONE = 0.d0
    ATWO = 0.d0
    do k=1,3
        AONE(k,k) = 2.d0
        ATWO(k,1:3) = 1.d0
    end do
    do k=4,6
        AONE(k,k) = 1.d0
    end do

    ! Initialize required matrices
    DDSDDE = 0.d0

    ! N-R Loop parameters
    ERROR = Y
    TOL=(10.d0**(-6))*Y
    MAXIT=30

    ! Deviatoric strain increment and deviatoric stress.
    DDEVST(1:6)= DSTRAN(1:6)-(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)/3.d0
    DEVS = STRESS-(STRESS(1)+STRESS(2)+STRESS(3))*(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)/3.d0

    ! Elastic predictors
    DEVSS = DEVS + ( E/(1+xnu) )*DDEVST
    SES = sqrt( 1.5d0* ( (DEVSS(1)**2) + (DEVSS(2)**2) + (DEVSS(3)**2) + &
        (DEVSS(4)**2) + (DEVSS(5)**2) + (DEVSS(6)**2) ))

    ! Newton-Raphson Loop. Solve for plastic strain increment.
    DEP = 0.d0
    if (SES*EDOT == 0)then
        DEP = 0.d0
    else
        do while (ERROR > TOL .and. NIT < MAXIT)
            T = Y*((1 +((EP+DEP)/xe))**(1/n))*((DEP/(DTIME*EDOT))**(1/m))
            F = SES - T - 1.5d0*(E/(1+xnu))*DEP
            DFDE = -1.5d0*(E/(1+xnu)) - T/(n*(xe+EP+DEP)) -T/(m*DEP)
            DEP = DEP - F/DFDE

            if (DEP<0) then
                DEP = DEP/10.d0
            end if

            ERROR = abs(F)
            NIT = NIT + 1
        end do
    end if

    ! Update Stress
    if (SES>0) then
        BETA = 1-1.5d0*E*DEP/((1+xnu)*SES)
    else
        BETA = 1.d0
    end if
    STRESS = BETA*DEVSS +&
    ((STRESS(1)+STRESS(2)+STRESS(3)) + &
    E*(DSTRAN(1)+DSTRAN(2)+DSTRAN(3))/(1-2.d0*xnu))*(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)/3.d0

    ! Tangent Stiffness
    if (SES*DEP>0) then
        BETA = 1-1.5d0*E*DEP/((1+xnu)*SES)
        GAMMAA = (1-BETA) + BETA*((1/(n*(xe+EP+DEP))) + (1/(m*DEP)))
        FACTOR = (E/(1+xnu))*(2.25*E*(DEP-1/GAMMAA)/((1+xnu)*(SES**3)))
    else
        BETA = 1.d0
        FACTOR = 0.d0
    end if

    ! Tangent stiffness
    ! Find Sij*Skl
    ATHREE = 0.d0
    do i=1,6
        do j=1,6
            ATHREE(i,j)=DEVSS(i)*DEVSS(j)
        end do
    end do
    DDSDDE = (E/(1+xnu))*BETA*(0.5d0*AONE-(1/3.d0)*ATWO) + FACTOR*ATHREE + &
    (E/(3.d0*(1-2.d0*xnu)))*ATWO

    ! Update state variables
    USTATV(1) = EP + DEP

        return
    end subroutine usermat_viscoplastic


