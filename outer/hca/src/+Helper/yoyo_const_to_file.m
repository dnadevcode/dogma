function [] = yoyo_const_to_file(yoyoBindingProb,ligandLength,idd)

 if nargin < 3
     idd = num2str(randi(1000));
 end
    [sortedSubseq, sortv, countATs, orderSeq,sortord,allDna] = sorted_NT(ligandLength);

    fd = fopen([idd,num2str(ligandLength),'_model.txt'],'w');
    for i=1:length(allDna)
        fprintf(fd, [int2nt(allDna{i}),' ',num2str(yoyoBindingProb(orderSeq(i))),'\n']);
    end

    fclose(fd);

end

