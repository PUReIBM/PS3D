!-------
! Independent module
!-------
! Some string manipulation functions
!-------
! Code:   Steinli
! Author: Chidambaram Narayanan
!         Nuclear Engineering Laboratory
!         ETH Zurich
!-------
MODULE string_funcs
  USE precision

  Implicit none
  Private
  Public:: newfilename, upper_to_lower, lower_to_upper &
          ,positionnml, number_of_lines                &
          ,append_integer_fixed, to_string

  INTERFACE number_of_lines
    MODULE PROCEDURE number_of_lines_unitno
    MODULE PROCEDURE number_of_lines_filename
  END INTERFACE

!-------
 Contains
!-------

  !----------
  !----------
  ! To count the number of lines in a file given the UNIT number
  ! of the file.
  !----------
  !----------
  FUNCTION number_of_lines_unitno(unitno)
    Integer:: number_of_lines_unitno
    Integer:: unitno

    Integer:: count
    Character:: c

    count = 0
    Do
      Read(unitno,"(A)",END=101)c
      count = count + 1
    Enddo
    101 Continue

    number_of_lines_unitno = count
  END FUNCTION number_of_lines_unitno


  !----------
  !----------
  ! To count the number of lines in a file given the name
  ! of the file.
  !----------
  !----------
  FUNCTION number_of_lines_filename(filename)
    Integer:: number_of_lines_filename
    Character(LEN=*):: filename

    Integer:: count, unitno
    Character:: c

    unitno = 102
    open(unit=unitno,file=filename,form="formatted")

    count = 0
    Do
      Read(unitno,"(A)",END=101)c
      count = count + 1
    Enddo
    101 Continue

    number_of_lines_filename = count
  END FUNCTION number_of_lines_filename

  !----------
  !----------
  ! To position the file pointer to read a particular NAMELIST.
  ! Can be used to make an input file such that the order of
  ! the NAMELISTS can be arbitrary.
  !----------
  !----------
  FUNCTION positionnml(unitno,nmlname)
    Logical:: positionnml
    Integer, Intent(in):: unitno
    Character(LEN=*), Intent(in):: nmlname

    Integer:: lengthnml, amppos
    Character(LEN=80):: line
    Character(LEN=LEN_TRIM(nmlname)+1):: nmlampstnd

    nmlampstnd = "&"//nmlname
    nmlampstnd = upper_to_lower(nmlampstnd)
    lengthnml = LEN(nmlampstnd)

    Rewind(unit=unitno)
    DO
      Read(unitno,'(A80)',END=999)line
      line = upper_to_lower(line)
      amppos = INDEX(line,nmlampstnd)
      If (amppos.ne.0) then
        BACKSPACE(unitno)
        positionnml = .true.
        Exit
      Endif
    ENDDO
    Return

    !-------
    ! the requested NAMELIST does not exist
    !-------
    999 positionnml = .false.
  END FUNCTION positionnml

  !----------
  !----------
  ! To append an integer tag to a string template
  !   - preceeding zeros in the intergers are removed
  !   - the second integer is optional
  !   - a hyphen '-' is inserted between the two intergers
  !   - returns a string of length 80 => TRIM command is
  !     recommended by the calling routine.
  !----------
  !----------
  FUNCTION newfilename(template,itnb,itne)
    Character(LEN=80):: newfilename
    Character(LEN=*), Intent(in):: template
    Integer, Intent(in):: itnb
    Integer, Intent(in), Optional:: itne

    Integer:: i, j, tmp, ic, ipos
    Integer:: sizeb, sizee, zero, tlen
    Integer, Dimension(:), Allocatable:: strb, stre

    newfilename(:) = " "

    If (itnb.ne.0) then
      sizeb = LOG10(REAL(itnb)) + 1
    Else
      sizeb = 1
    Endif
      
    allocate(strb(sizeb))

    !-------
    ! split integer into array of integers
    !-------
    tmp = itnb
    Do i = 1,sizeb
      j = sizeb - i + 1
      ic = MOD(tmp,10)
      strb(j) = ic
      tmp = tmp/10
    Enddo

    !-------
    ! create new file name
    !-------
    tlen = LEN_TRIM(template)

    Do ipos = 1, tlen
      newfilename(ipos:ipos) = template(ipos:ipos)
    Enddo

    !-------
    ! append the first iteration count
    !-------
    zero = IACHAR("0")

    Do i = 1, sizeb
      ipos = tlen + i
      tmp = zero + strb(i)
      newfilename(ipos:ipos) = ACHAR(tmp)
    Enddo

    !-------
    ! If ending integer is present
    !-------
    IF (PRESENT(itne)) then
      If (itne.ne.0) then
        sizee = LOG10(REAL(itne)) + 1
      Else
        sizee = 1
      Endif

      allocate(stre(sizee))

      !-------
      ! split integer into array of integers
      !-------
      tmp = itne
      Do i = 1,sizee
        j = sizee - i + 1
        ic = MOD(tmp,10)
        stre(j) = ic
        tmp = tmp/10
      Enddo

      !-------
      ! append a dash
      !-------
      ipos = ipos+1
      newfilename(ipos:ipos) = "-"

      !-------
      ! append the second iteration count
      !-------
      Do i = 1, sizee
        ipos = ipos + 1
        tmp = zero + stre(i)
        newfilename(ipos:ipos) = ACHAR(tmp)
      Enddo
    ENDIF
  END FUNCTION newfilename

  !----------
  !----------
  ! To append an integer tag to a string template
  !   - fixed string size given by ndigits.
  !   - returns a string of length 80 => TRIM command is
  !     recommended by the calling routine.
  !----------
  !----------
  FUNCTION append_integer_fixed(template,intval,ndigits)
    Character(LEN=80):: append_integer_fixed
    Character(LEN=*), Intent(in):: template
    Integer, Intent(in):: intval, ndigits

    Integer:: i, j, tmp, ic, ipos, size_int
    Integer:: zero, tlen
    Integer, Dimension(ndigits):: isuffix
    Integer, Dimension(:), Allocatable:: str_int

    append_integer_fixed(:) = " "

    If (intval.ne.0) then
      size_int = LOG10(REAL(intval)) + 1
    Else
      size_int = 1
    Endif

    Allocate(str_int(size_int))

    !-------
    ! split integer into array of integers
    !-------
    tmp = intval
    Do i = 1,size_int
      j = size_int - i + 1
      ic = MOD(tmp,10)
      str_int(j) = ic
      tmp = tmp/10
    Enddo

    isuffix(:) = 0
    Do i = 1,size_int
      isuffix(ndigits-size_int+i) = str_int(i)
    Enddo

    !-------
    ! create new file name
    !-------
    tlen = LEN_TRIM(template)

    Do ipos = 1, tlen
      append_integer_fixed(ipos:ipos) = template(ipos:ipos)
    Enddo

    !-------
    ! append the integer
    !-------
    zero = IACHAR("0")

    Do i = 1, ndigits
      ipos = tlen + i
      tmp = zero + isuffix(i)
      append_integer_fixed(ipos:ipos) = ACHAR(tmp)
    Enddo
  END FUNCTION append_integer_fixed

  !----------
  !----------
  ! To convert a string from upper to lower case
  !   - leaves non-alphabet elements unchanged
  !   - returns a string of length 80 => TRIM command is
  !     recommended for the calling routine.
  !----------
  !----------
  FUNCTION upper_to_lower(string)
    Character(LEN=80):: upper_to_lower
    Character(LEN=*):: string

    Integer:: i, strlen
    Integer, parameter:: aAdiff = IACHAR("a")-IACHAR("A")
   
    strlen = LEN_TRIM(string)

    Do i = 1,strlen
      if ("A" <= string(i:i) .AND. string(i:i) <= "Z") then
        upper_to_lower(i:i) = ACHAR(IACHAR(string(i:i)) + aAdiff)
      else
        upper_to_lower(i:i) = string(i:i)
      endif
    Enddo
    upper_to_lower(strlen+1:) = " "
  END FUNCTION upper_to_lower

  !----------
  !----------
  ! To convert a string from lower to upper case
  !   - leaves non-alphabet elements unchanged
  !   - returns a string of length 80 => TRIM command is
  !     recommended for the calling routine.
  !   - January 14, 2002
  ! (** not tested **)
  !----------
  !----------
  FUNCTION lower_to_upper(string)
    Character(LEN=80):: lower_to_upper
    Character(LEN=*):: string

    Integer:: i, strlen
    Integer, parameter:: aAdiff = IACHAR("a")-IACHAR("A")
   
    strlen = LEN_TRIM(string)

    Do i = 1,strlen
      if ("a" <= string(i:i) .AND. string(i:i) <= "z") then
        lower_to_upper(i:i) = ACHAR(IACHAR(string(i:i)) - aAdiff)
      else
        lower_to_upper(i:i) = string(i:i)
      endif
    Enddo
    lower_to_upper(strlen+1:) = " "
  END FUNCTION lower_to_upper


SUBROUTINE to_string(number, string, size)
	IMPLICIT NONE
	INTEGER number
	CHARACTER *(*) string
	INTEGER size
	CHARACTER*100 temp
	INTEGER local
	INTEGER last_digit
	INTEGER i

	local = number
	i = 0

	DO
	    last_digit = MOD(local,10)
	    local = local/10
	    i = i + 1
	    temp(i:i) = CHAR(last_digit + ICHAR('0'))
		IF (local==0) EXIT
	ENDDO
	
	size = i

	do i = 1, size
	    string(size-i+1:size-i+1) = temp(i:i)
	ENDDO
END SUBROUTINE to_string

INTEGER FUNCTION string_len(string)
	CHARACTER*(*) string

	CHARACTER*1 space
	PARAMETER (space = ' ')
	INTEGER i

	i = LEN(string)

	DO
		IF ((string(i:i).eq.space).and.(i.gt.1)) THEN
			i = i - 1
		ELSE 
			EXIT
		ENDIF
	ENDDO

	IF ((i.eq.1).and.(string(i:i).eq.space)) THEN
		string_len = 0
	ELSE
		string_len = i
	ENDIF
END FUNCTION string_len


END MODULE string_funcs

!-------
!PROGRAM teststrings
!  Use string_funcs
!
!  Integer:: i1, i2
!  Character(LEN=20):: template
!
!  i1 = 10; i2 = 12
!  template="evals"
!  write(*,*)TRIM(newfilename(template,i1,i2))
!  write(*,*)TRIM(newfilename(template,i1))
!END PROGRAM teststrings
!-------

!-------
! List of functions
!-------
! o Function positionnml
! o Function newfilename
! o Function upper_to_lower
! o Function lower_to_upper
!-------
