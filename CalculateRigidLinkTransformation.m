function TL = CalculateRigidLinkTransformation(Master_Node, Slave_Node)

X1 = Master_Node.Xcoord; Y1 = Master_Node.Ycoord;
X2 = Slave_Node.Xcoord;  Y2 = Slave_Node.Ycoord;

deltaX = X2 - X1;
deltaY = Y2 - Y1;

TL = [ 1 0 -deltaY;
       0 1  deltaX;
       0 0  1    ];
   

