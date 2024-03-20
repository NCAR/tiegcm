module matutil_module
! Perform direct calculations of small matrices (rank <= 4).
! Copied from https://fortranwiki.org/fortran/show/Matrix+inversion

  implicit none

  contains
!-------------------------------------------------------------------
  pure function matdet2(A) result(d)
! Calculate the determinant of the matrix

    real,dimension(2,2),intent(in) :: A
    real :: d

    d = A(1,1)*A(2,2) - A(1,2)*A(2,1)

  end function matdet2
!-------------------------------------------------------------------
  pure function matdet3(A) result(d)
! Calculate the determinant of the matrix

    real,dimension(3,3),intent(in) :: A
    real :: d

    d = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
      - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
      + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

  end function matdet3
!-------------------------------------------------------------------
  pure function matdet4(A) result(d)
! Calculate the determinant of the matrix

    real,dimension(4,4),intent(in) :: A
    real :: d

    d = A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))) &
      - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))) &
      + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))) &
      - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

  end function matdet4
!-------------------------------------------------------------------
  pure function matadj2(A) result(B)
! Calculate the adjugate of the matrix

    real,dimension(2,2),intent(in) :: A
    real,dimension(2,2) :: B

    B(1,1) =  A(2,2)
    B(2,1) = -A(2,1)
    B(1,2) = -A(1,2)
    B(2,2) =  A(1,1)

  end function matadj2
!-------------------------------------------------------------------
  pure function matadj3(A) result(B)
! Calculate the adjugate of the matrix

    real,dimension(3,3),intent(in) :: A
    real,dimension(3,3) :: B

    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))

  end function matadj3
!-------------------------------------------------------------------
  pure function matadj4(A) result(B)
! Calculate the adjugate of the matrix

    real,dimension(4,4),intent(in) :: A
    real,dimension(4,4) :: B

    B(1,1) = A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
    B(2,1) = A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
    B(3,1) = A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    B(4,1) = A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    B(1,2) = A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
    B(2,2) = A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
    B(3,2) = A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    B(4,2) = A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    B(1,3) = A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
    B(2,3) = A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
    B(3,3) = A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
    B(4,3) = A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
    B(1,4) = A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
    B(2,4) = A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    B(3,4) = A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
    B(4,4) = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

  end function matadj4
!-------------------------------------------------------------------
  pure function matinv2(A) result(B)
! Calculate the inverse of the matrix

    real,dimension(2,2),intent(in) :: A
    real,dimension(2,2) :: B

    B = matadj2(A)/matdet2(A)

  end function matinv2
!-------------------------------------------------------------------
  pure function matinv3(A) result(B)
! Calculate the inverse of the matrix

    real,dimension(3,3),intent(in) :: A
    real,dimension(3,3) :: B

    B = matadj3(A)/matdet3(A)

  end function matinv3
!-------------------------------------------------------------------
  pure function matinv4(A) result(B)
! Calculate the inverse of the matrix

    real,dimension(4,4),intent(in) :: A
    real,dimension(4,4) :: B

    B = matadj4(A)/matdet4(A)

  end function matinv4
!-------------------------------------------------------------------
endmodule matutil_module
