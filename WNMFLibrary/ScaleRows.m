function G=ScaleRows(G)
  [m1]=size(G);
  
  for i=1:m1
    seg=sum(G(i,:));
    G(i,:)=G(i,:)/seg;
  end
end