#ifndef _sz3_utils_hpp
#define _sz3_utils_hpp

template<class T, class T_eb>
class VariableEBLinearQuantizer{
public:
    VariableEBLinearQuantizer(int r = 32768) : radius(r) {}

    int get_radius() const { return radius; }

    // quantize the data with a prediction value, and returns the quantization index
    int quantize(T data, T pred, T_eb eb) {
        if(eb == 0) return 0;
        T diff = data - pred;
        int quant_index = (int) (fabs(diff) / eb) + 1;
        if (quant_index < this->radius * 2) {
            quant_index >>= 1;
            int half_index = quant_index;
            quant_index <<= 1;
            int quant_index_shifted;
            if (diff < 0) {
                quant_index = -quant_index;
                quant_index_shifted = this->radius - half_index;
            } else {
                quant_index_shifted = this->radius + half_index;
            }
            T decompressed_data = pred + quant_index * eb;
            if (fabs(decompressed_data - data) > eb) {
                return 0;
            } else {
                return quant_index_shifted;
            }
        } else {
            return 0;
        }
    }

    // quantize the data with a prediction value, and returns the quantization index and the decompressed data
    int quantize_and_overwrite(T &data, T pred, T_eb eb) {
        if(eb == 0){
            unpred.push_back(data);
            return 0;
        }
        // if(fabs(data + 14.3927) < 0.0001){
        //     std::cout << data << " " << pred << " " << eb << std::endl;
        // }
        T diff = data - pred;
        int quant_index = (int) (fabs(diff) / eb) + 1;
        if (quant_index < this->radius * 2) {
            quant_index >>= 1;
            int half_index = quant_index;
            quant_index <<= 1;
            int quant_index_shifted;
            if (diff < 0) {
                quant_index = -quant_index;
                quant_index_shifted = this->radius - half_index;
            } else {
                quant_index_shifted = this->radius + half_index;
            }
            T decompressed_data = pred + quant_index * eb;
            // if(fabs(data + 14.3927) < 0.0001){
            //     std::cout << decompressed_data << ", err = " << decompressed_data - data << std::endl;
            // }
            if (fabs(decompressed_data - data) > eb) {
                unpred.push_back(data);
                return 0;
            } else {
                data = decompressed_data;
                return quant_index_shifted;
            }
        } else {
            unpred.push_back(data);
            return 0;
        }
    }

    int quantize_and_overwrite(T ori, T pred, T_eb eb, T &dest) {
        if(eb == 0){
            unpred.push_back(ori);
            dest = ori;
            return 0;
        }
        T diff = ori - pred;
        int quant_index = (int) (fabs(diff) / eb) + 1;
        if (quant_index < this->radius * 2) {
            quant_index >>= 1;
            int half_index = quant_index;
            quant_index <<= 1;
            int quant_index_shifted;
            if (diff < 0) {
                quant_index = -quant_index;
                quant_index_shifted = this->radius - half_index;
            } else {
                quant_index_shifted = this->radius + half_index;
            }
            T decompressed_data = pred + quant_index * eb;
            if (fabs(decompressed_data - ori) > eb) {
                unpred.push_back(ori);
                dest = ori;
                return 0;
            } else {
                dest = decompressed_data;
                return quant_index_shifted;
            }
        } else {
            unpred.push_back(ori);
            dest = ori;
            return 0;
        }
    }

    // recover the data using the quantization index
    T recover(T pred, int quant_index, T_eb eb) {
        if (quant_index) {
            return recover_pred(pred, quant_index, eb);
        } else {
            return recover_unpred();
        }
    }

    T recover_pred(T pred, int quant_index, T_eb eb) {
        return pred + 2 * (quant_index - this->radius) * eb;
    }

    T recover_unpred() {
        return unpred[index++];
    }

    // required function in Quantizer interface
    int quantize(T data, T pred) {
        return 0;
    }

    // required function in Quantizer interface
    int quantize_and_overwrite(T &data, T pred) {
        return 0;
    }

    // required function in Quantizer interface
    T recover(T pred, int quant_index){
        return 0;
    }

    size_t size_est() {
        return unpred.size() * sizeof(T);
    }

    void save(unsigned char *&c) const {
        // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
        c[0] = 0b00000010;
        c += 1;
        // std::cout << "saving eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
        *reinterpret_cast<int *>(c) = this->radius;
        c += sizeof(int);
        *reinterpret_cast<size_t *>(c) = unpred.size();
        c += sizeof(size_t);
        memcpy(c, unpred.data(), unpred.size() * sizeof(T));
        c += unpred.size() * sizeof(T);
        // std::cout << "unpred size = " << unpred.size() << std::endl;
    };

    void load(const unsigned char *&c, size_t &remaining_length) {
        assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
        c += sizeof(uint8_t);
        this->radius = *reinterpret_cast<const int *>(c);
        c += sizeof(int);
        size_t unpred_size = *reinterpret_cast<const size_t *>(c);
        c += sizeof(size_t);
        this->unpred = std::vector<T>(reinterpret_cast<const T *>(c), reinterpret_cast<const T *>(c) + unpred_size);
        c += unpred_size * sizeof(T);
        // reset index
        index = 0;
    }

    void clear() {
        // std::cout << "unpred size = " << unpred.size() << std::endl;
        unpred.clear();
        index = 0;
    }

    virtual void postcompress_data() {}

    virtual void postdecompress_data() {}

    virtual void precompress_data() {}

    virtual void predecompress_data() {}


private:
    std::vector<T> unpred;
    size_t index = 0; // used in decompression only
    int radius; // quantization interval radius
};

template<class T>
inline T interp_linear(T a, T b) {
    return (a + b) / 2;
}

template<class T>
inline T interp_quad_1(T a, T b, T c) {
    return (3 * a + 6 * b - c) / 8;
}

template<class T>
inline T interp_quad_2(T a, T b, T c) {
    return (-a + 6 * b + 3 * c) / 8;
}

template<class T>
inline T interp_cubic(T a, T b, T c, T d) {
    return (-a + 9 * b + 9 * c - d) / 16;
}

#endif
