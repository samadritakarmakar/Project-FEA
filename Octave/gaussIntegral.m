function I= gaussIntegral(ch,degree,const,coordinates, stableCh)
    

    file=strcat('Data2/n',num2str(degree));
    load (file,'-mat');
    

    for i=1:degree
     for j=1:degree   
        I1=data(i,2)*data(j,2)*F(data(i,3),data(j,3),const, coordinates, ch, stableCh);
        if (i==1 && j==1)
            I =zeros(size(I1));
        end
        I=I+I1;
     end
    end
end
