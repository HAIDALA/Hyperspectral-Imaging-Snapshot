function rank = GetCaveImageRank(ImageName,ImagesArray,RankArray)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

j=find(ImagesArray==ImageName);
rank=RankArray(j);

end

