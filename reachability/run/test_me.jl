function test_main()
    print("Hello World\n")
end

# allows this script to be compiled into and executable
Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    test_main()  # call your program's logic.
    return 0
end
