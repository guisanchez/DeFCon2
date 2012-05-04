! Copyright (c) 2012 Joseph A. Levin
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or 
! substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.


! FSON MODULE 
!
! File:   fson.f95
! Author: Joseph A. Levin
!
! Created on March 10, 2012, 5:24 AM
!

      program example1

    ! Typical usage should only require an explicit use of the fson module.
    ! The other modules will be used privatley by fson as required.  
      use fson
      real a
      integer i
!     declare a pointer variable.  Always use a pointer with fson_value.
      type(fson_value), pointer :: value

    ! parse the json file
      value => fson_parse("input.js")

    ! print the parsed data to the console
C    call fson_print(value)

      call fson_get(value, "time.end", a)
      write(*,*) 'Duración de la simulación ',a
      
      call fson_get(value, "params.tanphi", a)
      write(*,*) 'Ángulo de equilibrio=', atan(a), 'rad. (' ,
     +     atan(a)*180.0/3.14,'º)'    
      
      call fson_get(value, "preproc_mesh",i)
      write(*,*) 'Preprocesado de la malla (1=Si/0=No)', i

C      call fson_get(value, "phoneNumber.[1].number", achar)
C      write(*,*) 'Y su telefono es ', achar

    ! extract data from the parsed value        

    ! clean up 
      call fson_destroy(value)



      end program example1
