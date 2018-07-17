# Sum_column.awk -
# Example: awk -f sum2.awk COL=3 < numbers.txt 
{ sum+= $COL}
END { print sum }