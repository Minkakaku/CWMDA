find . -type f -path "./*/outs/filtered_feature_bc_matrix/*" \
  -exec cp --parents {} /你的/分析文件夹/ \;