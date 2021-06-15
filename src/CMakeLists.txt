add_library (cpsz SHARED 
	sz_compression_utils.cpp
	sz_compress_cp_preserve_2d.cpp
	sz_compress_cp_preserve_3d.cpp
	sz_compress_3d.cpp
	sz_compress_pwr.cpp
	sz_decompression_utils.cpp
	sz_decompress_cp_preserve_2d.cpp
	sz_decompress_cp_preserve_3d.cpp
	sz_decompress_3d.cpp
	sz_decompress_pwr.cpp
	sz_huffman.cpp
	sz_lossless.cpp)
target_include_directories(cpsz PRIVATE ${ZSTD_INCLUDES})
target_include_directories(cpsz PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(cpsz ${ZSTD_LIBS})
file(GLOB header_files ${PROJECT_SOURCE_DIR}/include/*.hpp)

install(FILES ${header_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
install (TARGETS cpsz
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)
