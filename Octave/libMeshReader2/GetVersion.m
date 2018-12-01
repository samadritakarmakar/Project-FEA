function [i Version]=GetVersion(i, FileContents, WordStart, WordEnd)
    i=i+1
    Version=GetWord(i, FileContents, WordStart, WordEnd);
end
