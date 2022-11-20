function G=ScaleColumns(G)
  [m1,m2]=size(G);
  
  for i=1:m2
    seg=sum(G(:,i));
    G(:,i)=G(:,i)/seg;
  end
end