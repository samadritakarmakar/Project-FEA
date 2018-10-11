function pos=SkipToNextLine(FileContents, pos, FileLength)
    do
        strng = FileContents(pos);
        if (~strcmp(strng,sprintf('\n')))
            pos=pos+1;
        end
    until(FileContents(pos)==sprintf('\n') || pos==FileLength)
end
 
