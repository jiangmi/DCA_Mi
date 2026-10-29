# test/unit/phys/dca_step/cluster_solver

dca_add_gtest(self_energy_mixer_test
        FAST
        GTEST_MAIN
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS     ${DCA_LIBS}
        )

add_subdirectory(ctaux)
add_subdirectory(ctint)
add_subdirectory(shared_tools)
add_subdirectory(thread_qmci)
