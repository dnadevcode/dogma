function [ probBinding] = cb_theory(ntIntSeq, NETROPSINconc,YOYO1conc,yoyo1BindingConstant,values, untrustedRegion,isnetrop)
    % cb_theory
    % Computing cb theory
    %     Args:
    %         ntSeq, NETROPSINconc,YOYO1conc,yoyo1BindingConstant,values, untrustedRegion,isnetrop
    % 
    %     Returns:
    %         probBinding: probability binding vector
    %

    if nargin < 7
        isnetrop = 1;
    end

    % alternative theory, should we keep this?
    % if nargin < 6 % if no method was selected // Need to change this! will produce error now..
    %     NETROPSINconc = 6E-6;
    %     YOYO1conc = 4E-8;
    %     yoyo1BindingConstant =  [1E8];
    %     %netropsinBindingConstant = [ 5E5 1E8 ];
    %     untrustedRegion = 1000;
    %     W = 100.0;
    %     S = 8.0;
    %     netropsinBindingConstant = [S^4 W*S^3 (S^2)*(W^2) S*W^3 W^4];
    % end;
    
        

    % so this depends on ovrelap length.
    if untrustedRegion > 0
        ntIntSeq = [ntIntSeq(end-untrustedRegion+1:end) ntIntSeq ntIntSeq(1:untrustedRegion)];
    end

    nSize = size(ntIntSeq,2);

    % This will be our output vector
    yoyoBinding = zeros(nSize,1);
    netropsinBinding = zeros(nSize,1);

    


    % First, based on sequence ntIntSeq, we want to calculate which binding are
    % possible for each ligand. YOYO1 can bind always, while Netrospin is more
    % probable to bind if there are G's or C's
    ntNetrospinSeq = zeros(1,nSize);
    for i=nSize-4:-1:1 %1 = 'A' %2 = 'C', 3 = 'G', 4 = 'T'
        ntNetrospinSeq(i) =values(ntIntSeq(i),ntIntSeq(i+1),ntIntSeq(i+2),ntIntSeq(i+3));
    end


    % First, based on sequence ntIntSeq, we want to calculate which binding are
    % possible for each ligand. YOYO1 can bind always, while Netrospin is more
    % % probable to bind if there are G's or C's
    % ntNetrospinSeq = zeros(1,nSize);
    % SW = zeros(1,nSize);
    % for i=nSize-4:-1:1 %1 = 'A' %2 = 'C', 3 = 'G', 4 = 'T'
    %     if ntIntSeq(i) == 1 || ntIntSeq(i) == 4
    %         SW(i) = 1;
    %     end   
    %     ntNetrospinSeq(i) = ntNetrospinSeq(i+1)-SW(i+4)+SW(i);
    % end

    % transfer matrix
    % transferMat =         [1 0 0 0 1 0 0 0 1;
    %                        1 0 0 0 1 0 0 0 1;
    %                        0 1 0 0 0 0 0 0 0;
    %                        0 0 1 0 0 0 0 0 0;
    %                        0 0 0 NETROPSINconc 0 0 0 0 0;
    %                        1 0 0 0 1 0 0 0 1;
    %                        0 0 0 0 0 1 0 0 0;
    %                        0 0 0 0 0 0 1 0 0;
    %                        0 0 0 0 0 0 0 YOYO1conc*yoyo1BindingConstant 0];
    %                                   
    % two choices for Netropsin, based on if there are G and C letters in a
    % given quadromer or not

    %choice = [transferMat(5,4)*netropsinBindingConstant(1) transferMat(5,4)*netropsinBindingConstant(2)];
    ntNetrospinSeq = NETROPSINconc.*ntNetrospinSeq;

    leftVec = zeros(nSize+1,9);
    rightVec = zeros(9,nSize+1);
    maxEltLeft = zeros(1, nSize);
    maxEltRight = zeros(1, nSize);

    % the initial state. Before this state, three states are allowed: there was
    % nothing before, the last element of netropsin ligand was bound, the last
    % element of YOYO1 ligand was bound
    leftVec(1,:) = [1 0 0 0 1 0 0 0 1];
    leftVec(1,:) = leftVec(1,:)./norm(leftVec(1,:));

    rightVec(:,nSize+1) = transpose([1 0 0 0 0 0 0 0 0]);
    %rightVec(:,nSize+1)  = rightVec(:,nSize+1) ./norm( rightVec(:,nSize+1) );

    yoyoConst = YOYO1conc*yoyo1BindingConstant;


    for i=1:nSize
         %transferMat(5,4)= choice(ntNetrospinSeq(i)); % change only netropsin binding prob
         %leftVec(i+1,:) = leftVec(i,:)*transferMat;
         leftVecPrev = leftVec(i,:);
         leftVecNext = [leftVecPrev(1)+leftVecPrev(2)+leftVecPrev(6), ...
                    leftVecPrev(3),  leftVecPrev(4),leftVecPrev(5)*ntNetrospinSeq(i),leftVecPrev(1)+leftVecPrev(2)+leftVecPrev(6),...
                    leftVecPrev(7),leftVecPrev(8),  leftVecPrev(9)*yoyoConst,leftVecPrev(1)+leftVecPrev(2)+leftVecPrev(6)];

         %transferMat(5,4)= choice(ntNetrospinSeq(nSize-i+1)) ;

         %rightVec(:,nSize+1-i) = transferMat*rightVec(:,nSize-i+2);

         rightVecPrev = rightVec(:,nSize-i+2);
         rightVecNext = [rightVecPrev(1)+rightVecPrev(5)+rightVecPrev(9);
                         rightVecPrev(1)+rightVecPrev(5)+rightVecPrev(9);
                         rightVecPrev(2); rightVecPrev(3);
                         rightVecPrev(4)*ntNetrospinSeq(nSize-i+1);
                         rightVecPrev(1)+rightVecPrev(5)+rightVecPrev(9);
                         rightVecPrev(6); rightVecPrev(7); rightVecPrev(8)*yoyoConst];

         maxEltLeft(i) = norm(leftVecNext);
         maxEltRight(nSize+1-i) = norm(rightVecNext);

         %leftVec(i+1,:) = leftVec(i+1,:)./maxEltLeft(i);
         leftVec(i+1,:) = leftVecNext./maxEltLeft(i);
         %rightVec(:,nSize+1-i) = rightVec(:,nSize+1-i)./maxEltRight(nSize+1-i);
         rightVec(:,nSize+1-i) = rightVecNext./maxEltRight(nSize+1-i);

    end

    maxVecDiv =  zeros(1,nSize);
    maxVecDiv(1) = maxEltLeft(1)/maxEltRight(1);

    for i=2:nSize
        maxVecDiv(i) = maxVecDiv(i-1)*maxEltLeft(i)/maxEltRight(i);
    end

    denominator = leftVec(1,:)*rightVec(:,1);

    oMat = diag([0,0,0,0,0,1,1,1,1]); % this selects yoyo1 probability binding vec.
    yoyoBinding(1) = leftVec(1,:)*oMat*rightVec(:,1)/denominator;

    for i=2:nSize
        yoyoBinding(i) = leftVec(i,:)*oMat*rightVec(:,i)*maxVecDiv(i-1)/denominator;
    end

    if untrustedRegion > 0
        yoyoBinding = yoyoBinding(untrustedRegion+1:end-untrustedRegion);
    end
    
    if isnetrop == 2
        oMat = diag([0,1,1,1,1,0,0,0,0]); % this selects netropsin probability binding vec.
        netropsinBinding(1) = leftVec(1,:)*oMat*rightVec(:,1)/denominator;

        for i=2:nSize
            netropsinBinding(i) = leftVec(i,:)*oMat*rightVec(:,i)*maxVecDiv(i-1)/denominator;
        end

        if untrustedRegion > 0
            netropsinBinding = netropsinBinding(untrustedRegion+1:end-untrustedRegion);
        end
        probBinding = [yoyoBinding netropsinBinding];
    else
        probBinding = yoyoBinding;
    end

end


