# add_executable (sz_test sz_compress_test.cpp)
# target_include_directories (sz_test PRIVATE ${ZSTD_INCLUDES})
# target_include_directories (sz_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
# target_link_libraries (sz_test cpsz)

# add_executable (sz_test_with_eb sz_compress_test_with_eb.cpp)
# target_include_directories (sz_test_with_eb PRIVATE ${ZSTD_INCLUDES})
# target_include_directories (sz_test_with_eb PRIVATE ${PROJECT_SOURCE_DIR}/include)
# target_link_libraries (sz_test_with_eb cpsz)

# add_executable (sz_compress_cp_preserve_2d_test sz_compress_cp_preserve_2d_test.cpp)
# target_include_directories (sz_compress_cp_preserve_2d_test PRIVATE ${ZSTD_INCLUDES})
# target_include_directories (sz_compress_cp_preserve_2d_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
# target_link_libraries (sz_compress_cp_preserve_2d_test cpsz)

# add_executable (sz_compress_cp_preserve_3d_test sz_compress_cp_preserve_3d_test.cpp)
# target_include_directories (sz_compress_cp_preserve_3d_test PRIVATE ${ZSTD_INCLUDES})
# target_include_directories (sz_compress_cp_preserve_3d_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
# target_link_libraries (sz_compress_cp_preserve_3d_test cpsz)

# add_executable (sz_compress_cp_preserve_3d_unstructured_test sz_compress_cp_preserve_3d_unstructured_test.cpp)
# target_include_directories (sz_compress_cp_preserve_3d_unstructured_test PRIVATE ${ZSTD_INCLUDES})
# target_include_directories (sz_compress_cp_preserve_3d_unstructured_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
# target_link_libraries (sz_compress_cp_preserve_3d_unstructured_test cpsz)

add_executable (sz_compress_1d_with_eb sz_compress_1d_with_eb.cpp)
target_include_directories (sz_compress_1d_with_eb PRIVATE ${ZSTD_INCLUDES})
target_include_directories (sz_compress_1d_with_eb PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (sz_compress_1d_with_eb cpsz)

add_executable (sz_compress_2d_with_eb sz_compress_2d_with_eb.cpp)
target_include_directories (sz_compress_2d_with_eb PRIVATE ${ZSTD_INCLUDES})
target_include_directories (sz_compress_2d_with_eb PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (sz_compress_2d_with_eb cpsz)

add_executable (sz_compress_3d_with_eb sz_compress_3d_with_eb.cpp)
target_include_directories (sz_compress_3d_with_eb PRIVATE ${ZSTD_INCLUDES})
target_include_directories (sz_compress_3d_with_eb PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (sz_compress_3d_with_eb cpsz)

# install (TARGETS sz_compress_cp_preserve_2d_test sz_compress_cp_preserve_3d_test sz_compress_cp_preserve_3d_unstructured_test sz_compress_2d_with_eb
#   RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
# )
install (TARGETS sz_compress_2d_with_eb sz_compress_3d_with_eb
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)