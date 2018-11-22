function Word=GetWord(i, FileContents, WordStart, WordEnd)
    Word=FileContents(WordStart(i):WordEnd(i)); 
end
