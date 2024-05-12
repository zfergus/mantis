add_executable(mantis_debugger
        main.cpp
)

target_link_libraries(mantis_debugger PRIVATE mantis utils)
