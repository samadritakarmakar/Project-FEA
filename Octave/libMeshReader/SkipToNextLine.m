function pos=SkipToNextLine(FileContents, pos, FileLength)
   while(FileContents(pos)~=sprintf('\n') && pos~=FileLength)
    strng = FileContents(pos);
        if (~strcmp(strng,sprintf('\n')))
            pos=pos+1;
        end
    end
    %{
    do
        strng = FileContents(pos);
        if (~strcmp(strng,sprintf('\n')))
            pos=pos+1;
        end
    until(FileContents(pos)==sprintf('\n') || pos==FileLength)
%}
end

