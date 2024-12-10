subroutine SPECIFIC_CSV2ARRAY(numelems, csvString, array)

integer, intent(in) :: numelems
character*(*) :: csvString
SPECIFIC_TYPEKIND, intent(out) :: array(*)

integer :: strLength, i, startPos, endPos

if (numelems .lt. 1) return

strLength = len(trim(csvString))
startPos = 1
do i = 1, numelems-1
    endPos = scan(csvString(startPos:strLength),',') + startPos - 2
    array(i) = SPECIFIC_CONVERTER(csvString(startPos:endPos))
    startPos = endPos + 2
enddo
array(numelems) = SPECIFIC_CONVERTER(csvString(startPos:strLength))

end subroutine SPECIFIC_CSV2ARRAY
