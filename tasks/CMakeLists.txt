
# Create a main program using the library

set(EXECUTABLES run_pid_dcm12 run_pid_dcm3)

foreach(EXE ${EXECUTABLES})
    add_executable(${EXE} ${EXE}.cpp)
    target_link_libraries(${EXE} Pid ${ROOT_LIBRARIES})
    target_include_directories(${EXE} PUBLIC ${CMAKE_SOURCE_DIR}/src)
    install(TARGETS ${EXE} RUNTIME DESTINATION bin)
endforeach()


if (${PID_BUNDLED_AT})
    include(AnalysisTree)

    add_executable(fill_pid fill_pid.cpp)
    target_link_libraries(fill_pid Pid PidAT AnalysisTreeBase AnalysisTreeInfra ${ROOT_LIBRARIES})
    target_include_directories(fill_pid PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/at_interface ${PROJECT_INCLUDE_DIRECTORIES})

    install(TARGETS fill_pid RUNTIME DESTINATION bin)
endif ()


