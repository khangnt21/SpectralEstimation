
# Ước lượng Phổ - Spectral Estimatition
Bài báo cáo sử dụng phương pháp ước lượng phổ Blackman-Tukey, đồng thời áp dụng thêm các loại của sổ để phân tích và đánh giá. 
Sử dụng công cụ MATLAB, và các hàm mô phỏng các phương pháp đã được cung cấp trong môn học để viết bài báo cáo này.
Các loại của sổ sử dụng trong báo cáo: 

1. Cửa sổ chữ nhật
2. Cửa sổ Hanning
3. Cửa sổ Hamming
4. Cửa sổ Blackman
5. Cửa sổ Barlett

Phân tích ảnh hưởng của của sổ lên phương pháp phân tích phổ Blackman-Tuckey.

Khi ta xét một tín hiệu dừng bậc hai $X_c(t)$ có hàm tự tương quan $R^c_{XX}{f}$ tương ứng với phổ công suất $S^c_{XX}(f)$ có dải thông hữu hạn $B$. 
($S^c_{XX}(f)$ là kết quả mong muốn mà quá trình ước lượng muốn hội tụ đến).

Quá trình ước lượng sẽ bắt đầu với việc lấy mẫu tín hiệu $X_c(t)$ với một tần số $f_{\text{sampling}} > 2B$ (thoả mãn điều kiện Nyquist) 
để thu được tín hiệu $X_d(n)$. Sử dụng các phương pháp ước lượng ta sẽ tính được phổ $S^d_{XX}(f)$. Quá trình ước lượng mình sẽ cố gắng hội tụ kết quả đáp án xấp xỉ 
$S^d_{XX}(f)$ về đáp án mong muốn $S^c_{XX}(f)$.

Phương pháp Blackman-Tuckey sẽ tiến hành ước lượng hàm tự tương quan $R^z_{XX}(n)$ 
sau đó lấy biển đổi Fourier của $R^z_{XX}(n)$ để thu được phổ công suất $P^z_{XX}(n)$.