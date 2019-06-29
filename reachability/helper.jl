
function showtypetree(T, level=0)
    println("\t" ^ level, T)
    for t in subtypes(T)
        showtypetree(t, level+1)
   end
end

let
    counter = 0
    global function disp(msg)
        counter += 1
        print("[$counter]: $msg\n")
    end
end

