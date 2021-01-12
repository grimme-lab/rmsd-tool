! SPDX-Identifier: LGPL-3.0-or-later
module targ
   use tomlf
   implicit none
   private

   public :: arg_type, targ_type
   public :: get_arguments, new_argument_parser, get_options

   type :: arg_type
      character(len=:), allocatable :: arg
   end type arg_type

   type :: opt_type
      character(len=:), allocatable :: name
      integer :: require
   end type opt_type

   type :: targ_type
      integer :: nopt
      type(opt_type), allocatable :: opt(:)
   contains
      procedure :: add_option
   end type targ_type

   integer, parameter :: initial_size = 16

contains


subroutine get_arguments(args)
   type(arg_type), allocatable, intent(out) :: args(:)
   integer :: iarg, length
   character(len=:), allocatable :: tmp

   allocate(args(command_argument_count()))

   do iarg = 1, size(args)
      call get_command_argument(iarg, length=length)
      allocate(character(len=length) :: tmp)
      if (length > 0) then
         call get_command_argument(iarg, tmp)
      end if
      call move_alloc(tmp, args(iarg)%arg)
   end do

end subroutine get_arguments


subroutine get_options(self, args, table)
   type(targ_type), intent(in) :: self
   type(arg_type), allocatable, intent(inout) :: args(:)
   type(arg_type), allocatable :: tmp(:)
   type(toml_table), allocatable, intent(out) :: table
   type(toml_array), pointer :: array
   character(len=:), allocatable :: arg
   integer :: iarg, ipos, iopt, ii
   logical :: getopts
   table = toml_table()
   allocate(tmp(size(args)))
   getopts = .true.
   ipos = 0
   iarg = 0
   do while(iarg < size(args))
      iarg = iarg + 1
      call move_alloc(args(iarg)%arg, arg)
      if (getopts) then
         getopts = arg /= "--"
         if (.not.getopts) cycle
         call match_option(self, arg, iopt)
         if (iopt > 0) then
            select case(self%opt(iopt)%require)
            case(0)
               call set_value(table, self%opt(iopt)%name, .true.)
            case(1)
               iarg = iarg + 1
               if (iarg <= size(args)) then
                  call set_value(table, self%opt(iopt)%name, args(iarg)%arg)
               end if
            case default
               if (iarg + self%opt(iopt)%require <= size(args)) then
                  call add_array(table, self%opt(iopt)%name, array)
                  do ii = 1, self%opt(iopt)%require
                     call set_value(array, ii, args(iarg+ii)%arg)
                  end do
                  iarg = iarg + self%opt(iopt)%require
               end if
            end select
            cycle
         end if
      end if
      ipos = ipos + 1
      call move_alloc(arg, tmp(ipos)%arg)
   end do

   deallocate(args)
   allocate(args(ipos))
   do ipos = 1, size(args)
      call move_alloc(tmp(ipos)%arg, args(ipos)%arg)
   end do

end subroutine get_options


subroutine match_option(self, arg, iopt)
   type(targ_type), intent(in) :: self
   character(len=*), intent(in) :: arg
   integer, intent(out) :: iopt
   integer :: ii
   iopt = 0
   do ii = 1, self%nopt
      if (self%opt(ii)%name == arg(3:)) then
         iopt = ii
         exit
      end if
   end do
end subroutine match_option


function new_argument_parser(nopt) result(self)
   integer, intent(in), optional :: nopt
   type(targ_type) :: self

   self%nopt = 0
   if (present(nopt)) then
      allocate(self%opt(nopt))
   else
      allocate(self%opt(initial_size))
   end if

end function new_argument_parser


subroutine add_option(self, name, require)
   class(targ_type), intent(inout) :: self
   character(len=*), intent(in) :: name
   integer, intent(in), optional :: require
   integer :: m

   m = size(self%opt)
   if (self%nopt >= m) then
      call resize(self%opt, m + m/2 + 1)
   end if

   self%nopt = self%nopt + 1
   self%opt(self%nopt)%name = name
   if (present(require)) then
      self%opt(self%nopt)%require = require
   else
      self%opt(self%nopt)%require = 0
   end if

end subroutine add_option


subroutine resize(list, n)

   !> Array to be resized
   type(opt_type), allocatable, intent(inout), target :: list(:)

   !> New size of the list
   integer, intent(in) :: n

   type(opt_type), allocatable, target :: tmp(:)
   integer :: i


   if (allocated(list)) then
      call move_alloc(list, tmp)
      allocate(list(n))

      do i = 1, min(size(tmp), n)
         list(i) = tmp(i)
      end do

      deallocate(tmp)
   else
      allocate(list(n))
   end if

end subroutine resize


end module targ
