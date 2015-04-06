
botLeft = [-1, -.2, -.2]
size = [2.0, .4, .4]
n = [30, 12, 12] 
for i in range(n[0] + 1):
   for j in range(n[1] + 1):
      for k in range(n[2] + 1):
         x = botLeft[0] + float(i)/n[0]*size[0]
         y = botLeft[1] + float(j)/n[1]*size[1]
         z = botLeft[2] + float(k)/n[2]*size[2]
         print x, y, z
      #
   #
#
