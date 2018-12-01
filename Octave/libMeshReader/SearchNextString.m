function [strng pos]=SearchNextString(FileContents, pos, FileLength)
    i=1;
    while ((FileContents(pos)~=' ' && FileContents(pos)~=sprintf('\n') && pos~=FileLength) || i==1)
        
        CharCheck=FileContents(pos);
            if (CharCheck~=' ' || CharCheck~=sprintf('\n'))
                strng(i)=CharCheck;
                i=i+1;
            end  
        pos=pos+1;
    end
    %{
    do
        CharCheck=FileContents(pos);
        if (CharCheck~=' ' || CharCheck~=sprintf('\n'))
            strng(i)=CharCheck;
            i=i+1;
        end
        pos=pos+1;
    until(FileContents(pos)==' ' || FileContents(pos)==sprintf('\n') || pos==FileLength)
%}
end
