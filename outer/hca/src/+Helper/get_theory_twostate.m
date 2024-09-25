function [theorySeq, theoryStr] = get_theory_twostate(fastaName,gcSF,pxSize,nmpx,isC,sigma,kN,psf,cY, cN,kY,ligandLength)

    fasta = fastaread(fastaName);
    ntSeq = nt2int(fasta.Sequence);

    % cummulative sum of AT's. 
    numWsCumSum = cumsum((ntSeq == 1)  | (ntSeq == 4) );

    import CBT.SimpleTwoState.gen_simple_theory_px;
    [theorySeq] = gen_simple_theory_px(numWsCumSum,gcSF,pxSize,nmpx,isC,sigma,kN,psf,cY, cN,kY,ligandLength);

    if nargout >=2
        theoryStr{1}.rawBarcode = theorySeq;
        theoryStr{1}.rawBitmask =[];
        theoryStr{1}.length = length(theorySeq);
        theoryStr{1}.name = fasta.Header;
        theoryStr{1}.isLinearTF = ~isC;
    end
end

