#include "PairID.hpp"

PairID::PairID() : setBits(localBitVector) {}

PairID::PairID(const PairID& o) {

    const uintptr_t flag = o.setBits & flagMask;

    if (flag == ptrBitmap){

        Bitmap* setPtrBmp = new Bitmap;

        setPtrBmp->r = o.getConstPtrBitmap()->r;

        setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
    }
    else if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = o.getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        TinyBitmap t_bmp_cpy(t_bmp);

        t_bmp.detach();

        setBits = (reinterpret_cast<uintptr_t>(t_bmp_cpy.detach()) & pointerMask) | localTinyBitmap;
    }
    else setBits = o.setBits;
}

PairID::PairID(PairID&& o) : setBits(o.setBits) {

    o.setBits = localBitVector;
}

PairID::~PairID() {

    releaseMemory();
}

PairID& PairID::operator=(const PairID& o){

    if (this != &o) {

        const uintptr_t flag = o.setBits & flagMask;
        const uintptr_t this_flag = setBits & flagMask;

        if (flag == ptrBitmap){

            Bitmap* setPtrBmp = nullptr;

            if (this_flag == ptrBitmap) setPtrBmp = getPtrBitmap();
            else {

                releaseMemory();
                setPtrBmp = new Bitmap;
            }

            setPtrBmp->r = o.getConstPtrBitmap()->r;

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
        }
        else if (flag == localTinyBitmap){

            releaseMemory();

            uint16_t* setPtrTinyBmp = o.getPtrTinyBitmap();

            TinyBitmap t_bmp(&setPtrTinyBmp);
            TinyBitmap t_bmp_cpy(t_bmp);

            t_bmp.detach();

            setBits = (reinterpret_cast<uintptr_t>(t_bmp_cpy.detach()) & pointerMask) | localTinyBitmap;
        }
        else {

            releaseMemory();

            setBits = o.setBits;
        }
    }

    return *this;
}

PairID& PairID::operator=(PairID&& o){

    if (this != &o) {

        releaseMemory();

        setBits = o.setBits;
        o.setBits = localBitVector;
    }

    return *this;
}

bool PairID::operator==(const PairID& o) const {

    if (size() != o.size()) return false;

    PairID::const_iterator it(begin()), o_it(o.begin());
    const PairID::const_iterator it_end(end()), o_it_end(o.end());

    for (; (it != it_end) && (o_it != o_it_end); ++it, ++o_it){

        if (*it != *o_it) return false;
    }

    return ((it == it_end) && (o_it == o_it_end));
}

void PairID::clear(){

    releaseMemory();
    setBits = localBitVector;
}

size_t PairID::getSizeInBytes() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.getSizeInBytes() + sizeof(Bitmap) + sizeof(PairID);

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.getSizeInBytes();

        t_bmp.detach();

        return ret;
    }

    return sizeof(PairID); 
}

PairID PairID::operator|(const PairID& rhs) const {

    PairID lhs(*this);

    lhs |= rhs;
    return lhs;
}

PairID& PairID::operator|=(const PairID& rhs) {

    if (&rhs != this){

        if (cardinality() == 0) *this = rhs;
        else if (rhs.cardinality() != 0){
    
            const uintptr_t flag = setBits & flagMask;
            const uintptr_t r_flag = rhs.setBits & flagMask;

            if ((flag == r_flag) && ((flag == ptrBitmap) || (flag == localBitVector))) {

                if (flag == ptrBitmap) getPtrBitmap()->r |= rhs.getConstPtrBitmap()->r;
                else setBits |= rhs.setBits;
            }
            else {

                const size_t max_a = maximum(), min_a = minimum();
                const size_t max_b = rhs.maximum(), min_b = rhs.minimum();

                vector<uint32_t> new_ids;

                if ((min_a <= max_b) && (min_b <= max_a)) { // Check that range overlaps (both bitmaps must be non empty!)

                    const_iterator it = begin(min_b), it_end = end();
                    const_iterator r_it = rhs.begin(), r_it_end = rhs.end();

                    new_ids.reserve(cardinality() + rhs.cardinality());

                    while ((it != it_end) && (r_it != r_it_end)){

                        if (*it > *r_it){

                            new_ids.push_back(*r_it);
                            ++r_it;
                        }
                        else if (*it < *r_it) ++it;
                        else {

                            ++it;
                            ++r_it;
                        }
                    }

                    while (r_it != r_it_end){

                        new_ids.push_back(*r_it);
                        ++r_it;
                    }
                }
                else new_ids = rhs.toVector();

                addSortedVector(new_ids);
            }
        }
    }

    return *this;
}

PairID PairID::operator-(const PairID& rhs) const {

    PairID lhs(*this);

    lhs -= rhs;
    return lhs;
}

PairID& PairID::operator-=(const PairID& rhs) {

    const size_t a_card = cardinality();
    const size_t b_card = rhs.cardinality();

    if ((&rhs != this) && (a_card != 0) && (b_card != 0)) {

        const uintptr_t flag = setBits & flagMask;
        const uintptr_t r_flag = rhs.setBits & flagMask;

        if ((flag == r_flag) && ((flag == ptrBitmap) || (flag == localBitVector))) {

            if (flag == ptrBitmap) getPtrBitmap()->r -= rhs.getConstPtrBitmap()->r;
            if (flag == localBitVector) setBits -= (setBits & rhs.setBits & pointerMask);
        }
        else if (flag == localSingleInt) {

            if (rhs.contains(setBits >> shiftMaskBits)) clear();
        }
        else {

            const size_t max_a = maximum(), min_a = minimum();
            const size_t max_b = rhs.maximum(), min_b = rhs.minimum();

            if ((min_a <= max_b) && (min_b <= max_a)) { // Check that range overlaps (both bitmaps must be non empty!)

                const size_t log2_a = b_card * min(l_approximate_log2(a_card), 16UL);
                const size_t log2_b = a_card * min(l_approximate_log2(b_card), 16UL);

                const size_t min_ab_card = min(a_card + b_card, min(log2_a, log2_b));
                const size_t max_ab = max(min_a, min_b);

                const_iterator it = begin(max_ab), it_end = end();
                const_iterator r_it = rhs.begin(max_ab), r_it_end = rhs.end();

                vector<uint32_t> old_ids;

                old_ids.reserve(min(a_card, b_card));

                if (min_ab_card == (a_card + b_card)) {

                    while ((it != it_end) && (r_it != r_it_end)){

                        if (*it > *r_it) ++r_it;
                        else if (*it < *r_it) ++it;
                        else {

                            old_ids.push_back(*r_it);

                            ++it;
                            ++r_it;
                        }
                    }
                }
                else if (min_ab_card == log2_a){

                    while ((r_it != r_it_end) && (*r_it <= max_a)){

                        if (contains(*r_it)) old_ids.push_back(*r_it);

                        ++r_it;
                    }
                }
                else {

                    while ((it != it_end) && (*it <= max_b)){

                        if (rhs.contains(*it)) old_ids.push_back(*it);

                        ++it;
                    }
                }

                removeSortedVector(old_ids);
            }
        }
    }

    return *this;
}

PairID PairID::operator&(const PairID& rhs) const {

    PairID lhs(*this);

    lhs &= rhs;
    return lhs;
}

PairID& PairID::operator&=(const PairID& rhs) {
    
    if (&rhs != this){

        const size_t a_card = cardinality();
        const size_t b_card = rhs.cardinality();

        if ((a_card == 0) || (b_card == 0)) clear();
        else {

            const uintptr_t flag = setBits & flagMask;
            const uintptr_t r_flag = rhs.setBits & flagMask;

            if ((flag == r_flag) && ((flag == ptrBitmap) || (flag == localBitVector))) {

                if (flag == ptrBitmap) getPtrBitmap()->r &= rhs.getConstPtrBitmap()->r;
                if (flag == localBitVector) setBits &= rhs.setBits;
            }
            else {

                const size_t max_a = maximum(), min_a = minimum();
                const size_t max_b = rhs.maximum(), min_b = rhs.minimum();

                if ((min_a <= max_b) && (min_b <= max_a)) { // Check that range overlaps (both bitmaps must be non empty!)

                    const size_t log2_a = b_card * min(l_approximate_log2(a_card), 16UL);
                    const size_t log2_b = a_card * min(l_approximate_log2(b_card), 16UL);

                    const size_t min_ab_card = min(a_card + b_card, min(log2_a, log2_b));
                    const size_t max_ab = max(min_a, min_b);

                    PairID out_pid;

                    PairID::const_iterator a_it_s = begin(max_ab), a_it_e = end();
                    PairID::const_iterator b_it_s = rhs.begin(max_ab), b_it_e = rhs.end();

                    if (min_ab_card == (a_card + b_card)) {

                        while ((a_it_s != a_it_e) && (b_it_s != b_it_e)){

                            const uint32_t val_a = *a_it_s;
                            const uint32_t val_b = *b_it_s;

                            if (val_a == val_b){

                                out_pid.add(val_a);

                                ++a_it_s;
                                ++b_it_s;
                            }
                            else if (val_a < val_b) ++a_it_s;
                            else ++b_it_s;
                        }
                    }
                    else if (min_ab_card == log2_a){

                        while ((b_it_s != b_it_e) && (*b_it_s <= max_a)){

                            if (contains(*b_it_s)) out_pid.add(*b_it_s);

                            ++b_it_s;
                        }
                    }
                    else {

                        while ((a_it_s != a_it_e) && (*a_it_s <= max_b)){

                            if (rhs.contains(*a_it_s)) out_pid.add(*a_it_s);

                            ++a_it_s;
                        }
                    }

                    *this = move(out_pid);
                }
                else clear();
            }
        }
    }

    return *this;
}

size_t PairID::and_cardinality(const PairID& rhs) const {

    const size_t a_card = cardinality();
    const size_t b_card = rhs.cardinality();

    if ((a_card == 0) || (b_card == 0)) return 0;
    else {

        const uintptr_t flag = setBits & flagMask;
        const uintptr_t r_flag = rhs.setBits & flagMask;

        if ((flag == r_flag) && ((flag == ptrBitmap) || (flag == localBitVector))) {

            if (flag == ptrBitmap) return getPtrBitmap()->r.and_cardinality(rhs.getConstPtrBitmap()->r);
            
            return __builtin_popcountll(pointerMask & setBits & rhs.setBits);
        }
        else if (flag == localSingleInt) {

            return static_cast<size_t>(rhs.contains(setBits >> shiftMaskBits));
        }
        else if (r_flag == localSingleInt) {

            return static_cast<size_t>(contains(rhs.setBits >> shiftMaskBits));
        }
        else {

            const size_t max_a = maximum(), min_a = minimum();
            const size_t max_b = rhs.maximum(), min_b = rhs.minimum();

            if ((min_a <= max_b) && (min_b <= max_a)) { // Check that range overlaps (both bitmaps must be non empty!)

                const size_t log2_a = b_card * min(l_approximate_log2(a_card), 16UL);
                const size_t log2_b = a_card * min(l_approximate_log2(b_card), 16UL);

                const size_t min_ab_card = min(a_card + b_card, min(log2_a, log2_b));
                const size_t max_ab = max(min_a, min_b);

                size_t nb_shared = 0;

                PairID::const_iterator a_it_s = begin(max_ab), a_it_e = end();
                PairID::const_iterator b_it_s = rhs.begin(max_ab), b_it_e = rhs.end();

                if (min_ab_card == (a_card + b_card)) {

                    while ((a_it_s != a_it_e) && (b_it_s != b_it_e)){

                        const uint32_t val_a = *a_it_s;
                        const uint32_t val_b = *b_it_s;

                        if (val_a == val_b){

                            ++nb_shared;
                            ++a_it_s;
                            ++b_it_s;
                        }
                        else if (val_a < val_b) ++a_it_s;
                        else ++b_it_s;
                    }
                }
                else if (min_ab_card == log2_a){

                    while ((b_it_s != b_it_e) && (*b_it_s <= max_a)){

                        nb_shared += static_cast<size_t>(contains(*b_it_s));
                        ++b_it_s;
                    }
                }
                else {

                    while ((a_it_s != a_it_e) && (*a_it_s <= max_b)){

                        nb_shared += static_cast<size_t>(rhs.contains(*a_it_s));
                        ++a_it_s;
                    }
                }

                return nb_shared;
            }
            else return 0;
        }
    }
}

size_t PairID::and_cardinality(const PairID& rhs, const uint64_t min_shared) const {

    const size_t a_card = cardinality();
    const size_t b_card = rhs.cardinality();

    if ((a_card == 0) || (b_card == 0) || (min_shared == 0)) return 0;
    else {

        const uintptr_t flag = setBits & flagMask;
        const uintptr_t r_flag = rhs.setBits & flagMask;

        if ((flag == r_flag) && ((flag == ptrBitmap) || (flag == localBitVector))) {

            if (flag == ptrBitmap) return getPtrBitmap()->r.and_cardinality(rhs.getConstPtrBitmap()->r, min_shared);
            
            return __builtin_popcountll(pointerMask & setBits & rhs.setBits);
        }
        else if (flag == localSingleInt) return static_cast<size_t>(rhs.contains(setBits >> shiftMaskBits));
        else if (r_flag == localSingleInt) return static_cast<size_t>(contains(rhs.setBits >> shiftMaskBits));
        else {

            const size_t max_a = maximum(), min_a = minimum();
            const size_t max_b = rhs.maximum(), min_b = rhs.minimum();

            if ((min_a <= max_b) && (min_b <= max_a)) { // Check that range overlaps (both bitmaps must be non empty!)

                const size_t log2_a = b_card * min(l_approximate_log2(a_card), 16UL);
                const size_t log2_b = a_card * min(l_approximate_log2(b_card), 16UL);

                const size_t min_ab_card = min(a_card + b_card, min(log2_a, log2_b));
                const size_t max_ab = max(min_a, min_b);

                size_t nb_shared = 0;

                PairID::const_iterator a_it_s = begin(max_ab), a_it_e = end();
                PairID::const_iterator b_it_s = rhs.begin(max_ab), b_it_e = rhs.end();

                if (min_ab_card == (a_card + b_card)) {

                    while ((a_it_s != a_it_e) && (b_it_s != b_it_e) && (nb_shared < min_shared)){

                        const uint32_t val_a = *a_it_s;
                        const uint32_t val_b = *b_it_s;

                        if (val_a == val_b){

                            ++nb_shared;
                            ++a_it_s;
                            ++b_it_s;
                        }
                        else if (val_a < val_b) ++a_it_s;
                        else ++b_it_s;
                    }
                }
                else if (min_ab_card == log2_a){

                    while ((b_it_s != b_it_e) && (*b_it_s <= max_a) && (nb_shared < min_shared)){

                        nb_shared += static_cast<size_t>(contains(*b_it_s));
                        ++b_it_s;
                    }
                }
                else {

                    while ((a_it_s != a_it_e) && (*a_it_s <= max_b) && (nb_shared < min_shared)){

                        nb_shared += static_cast<size_t>(rhs.contains(*a_it_s));
                        ++a_it_s;
                    }
                }

                return nb_shared;
            }
            else return 0;
        }
    }
}

Roaring PairID::toRoaring() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getPtrBitmap()->r;
    else {

        Roaring b;

        const_iterator it = begin();
        const_iterator it_end = end();

        while (it != it_end){

            b.add(*it);
            ++it;
        }

        return b;
    }
}

vector<uint32_t> PairID::toVector() const {

    vector<uint32_t> v;

    v.reserve(cardinality());

    for (const uint32_t id : *this) v.push_back(id);

    return v;
}

void PairID::add(const size_t id) {

    uintptr_t flag = setBits & flagMask;

    if (flag == localSingleInt){

        const uintptr_t setBits_tmp = setBits >> shiftMaskBits;

        if (setBits_tmp != id){

            if ((setBits_tmp < maxBitVectorIDs) && (id < maxBitVectorIDs)){

                setBits = (1ULL << (setBits_tmp + shiftMaskBits)) | localBitVector;
            }
            else {

                TinyBitmap t_bmp;

                if (t_bmp.add(setBits_tmp)) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                else {

                    t_bmp.clear();

                    Bitmap* setPtrBmp = new Bitmap;

                    setPtrBmp->r.add(setBits_tmp);

                    setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
                }
            }
        

            flag = setBits & flagMask;
        }
    }

    if (flag == localBitVector){

        if (setBits == localBitVector) setBits = (id << shiftMaskBits) | localSingleInt;
        else if (id < maxBitVectorIDs) setBits |= 1ULL << (id + shiftMaskBits);
        else {

            uintptr_t setBits_tmp_tb = setBits >> shiftMaskBits;
            uintptr_t setBits_tmp_cr = setBits >> shiftMaskBits;

            TinyBitmap t_bmp;

            bool add_ok = true;

            for (size_t i = 0; (setBits_tmp_tb != 0) && add_ok; ++i, setBits_tmp_tb >>= 1) {

                if (setBits_tmp_tb & 0x1) add_ok = t_bmp.add(i);
            }

            if (add_ok) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            else {

                t_bmp.clear();

                Bitmap* setPtrBmp = new Bitmap;

                for (size_t i = 0; setBits_tmp_cr != 0; ++i, setBits_tmp_cr >>= 1) {

                    if (setBits_tmp_cr & 0x1) setPtrBmp->r.add(i);
                }

                setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            }
        }

        flag = setBits & flagMask;
    }

    if (flag == localTinyBitmap) {

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

        TinyBitmap t_bmp(&setPtrTinyBmp);

        if (t_bmp.add(id)) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        else {

            const size_t sz_t_bmp = t_bmp.size();

            size_t i = 0;

            uint32_t* values = new uint32_t[sz_t_bmp];

            Bitmap* setPtrBmp = new Bitmap;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) values[i++] = *it;

            t_bmp.clear();
            setPtrBmp->r.addMany(sz_t_bmp, values);

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) getPtrBitmap()->r.add(id); // flag == ptrBitmap
}

PairID PairID::forceRoaringInternal() const {

    PairID pid;
    Bitmap* setPtrBmp = new Bitmap;

    pid.setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
    setPtrBmp->r = toRoaring();

    return pid;
}

void PairID::addSortedVector(const vector<uint32_t>& v) { // Private, assumes vector is sorted

    if (v.empty()) return;

    uintptr_t flag = setBits & flagMask;

    size_t i = 0;

    const uint32_t v_back = v.back();
    const bool sameUpper16bits = ((v.front() >> 16) == (v_back >> 16));

    if (flag == localSingleInt){

        const uintptr_t setBits_tmp = setBits >> shiftMaskBits;

        if ((setBits_tmp < maxBitVectorIDs) && (v_back < maxBitVectorIDs)){

            setBits = (1ULL << (setBits_tmp + shiftMaskBits)) | localBitVector;
        }
        else {

            TinyBitmap t_bmp;

            if (sameUpper16bits && t_bmp.add(setBits_tmp)) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            else {

                t_bmp.clear();

                Bitmap* setPtrBmp = new Bitmap;

                setPtrBmp->r.add(setBits_tmp);

                setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            }
        }

        flag = setBits & flagMask;
    }

    if (flag == localBitVector){

        if ((setBits == localBitVector) && (v.size() == 1)) setBits = (v.front() << shiftMaskBits) | localSingleInt;
        else if (v_back < maxBitVectorIDs){

            for (const uint32_t id : v) setBits |= 1ULL << (id + shiftMaskBits);
        }
        else {

            uintptr_t setBits_tmp_tb = setBits >> shiftMaskBits;
            uintptr_t setBits_tmp_cr = setBits >> shiftMaskBits;

            TinyBitmap t_bmp;

            bool add_ok = sameUpper16bits;

            for (size_t j = 0; (setBits_tmp_tb != 0) && add_ok; ++j, setBits_tmp_tb >>= 1) {

                if (setBits_tmp_tb & 0x1) add_ok = t_bmp.add(j);
            }

            if (add_ok) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            else {

                Bitmap* setPtrBmp = new Bitmap;

                t_bmp.clear();

                for (size_t j = 0; setBits_tmp_cr != 0; ++j, setBits_tmp_cr >>= 1) {

                    if (setBits_tmp_cr & 0x1) setPtrBmp->r.add(j);
                }

                setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            }
        }

        flag = setBits & flagMask;
    }

    if (flag == localTinyBitmap) {

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

        TinyBitmap t_bmp(&setPtrTinyBmp);

        bool add_ok = sameUpper16bits && ((t_bmp.size() == 0) || ((v_back >> 16) == (*(t_bmp.begin()) >> 16)));

        for (i = 0; (i < v.size()) && add_ok; ++i) add_ok = t_bmp.add(v[i]);

        if (add_ok) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        else {

            Bitmap* setPtrBmp = new Bitmap;

            const size_t sz_t_bmp = t_bmp.size();

            if (sz_t_bmp != 0) {

                uint32_t* values = new uint32_t[sz_t_bmp];

                size_t j = 0;

                for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) values[j++] = *it;

                t_bmp.clear();

                setPtrBmp->r.addMany(sz_t_bmp, values);

                delete[] values;
            }

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            i -= static_cast<size_t>(i != 0);
        }
    }

    if (flag == ptrBitmap) {

        Bitmap* bitmap = getPtrBitmap();

        bitmap->r.addMany(v.size() - i, &v[i]);
        bitmap->r.runOptimize();
    }
}

void PairID::remove(const size_t pair_id) {

    uintptr_t flag = setBits & flagMask;

    if (flag == localBitVector){

        if (pair_id < maxBitVectorIDs) setBits &= ~(1ULL << (pair_id + shiftMaskBits));
    }
    else if (flag == localSingleInt){

        if (pair_id == (setBits >> shiftMaskBits)) setBits = localBitVector;
    }
    else if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

        TinyBitmap t_bmp(&setPtrTinyBmp);

        const bool rm_ok = t_bmp.remove(pair_id);

        if (rm_ok){

            const size_t card = t_bmp.size();

            if (card == 0){

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                clear();
            }
            else if (card == 1){

                const uint32_t l_pair_id = *(t_bmp.begin());

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;

                clear();
                add(l_pair_id);
            }
            else if (t_bmp.maximum() < maxBitVectorIDs){

                PairID new_uc;

                for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) new_uc.add(*it);

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;

                *this = move(new_uc);
            }
            else setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        }
        else {

            const size_t sz_t_bmp = t_bmp.size();

            uint32_t* values = new uint32_t[sz_t_bmp];

            size_t i = 0;

            Bitmap* setPtrBmp = new Bitmap;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) values[i++] = *it;

            t_bmp.clear();
            setPtrBmp->r.addMany(sz_t_bmp, values);

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) {

        Bitmap* bitmap = getPtrBitmap();

        bitmap->r.remove(pair_id);

        const size_t card = bitmap->r.cardinality();

        if (card == 0) clear();
        else if (card == 1){

            const uint32_t l_pair_id = bitmap->r.minimum();

            clear();
            add(l_pair_id);
        }
        else if (bitmap->r.maximum() < maxBitVectorIDs){

            PairID new_uc;

            const_iterator it = begin(), it_end = end();

            for (; it != it_end; ++it) new_uc.add(*it);

            *this = move(new_uc);
        }
        else if ((setBits & flagMask) == ptrBitmap) bitmap->r.runOptimize();
    }
}

void PairID::removeSortedVector(const vector<uint32_t>& v) {

    if (v.empty()) return;

    uintptr_t flag = setBits & flagMask;

    size_t i = 0;

    if (flag == localBitVector){

        uintptr_t mask = 0;

        for (size_t j = 0; j < v.size(); ++j){

            if (v[j] < maxBitVectorIDs) mask |= 1ULL << (v[j] + shiftMaskBits);
        }

        setBits &= ~mask;
    }
    else if (flag == localSingleInt){

        const uintptr_t id = setBits >> shiftMaskBits;

        if (id <= v[v.size() - 1]){

            for (size_t j = 0; j < v.size(); ++j){

                if (v[j] == id) setBits = localBitVector;
                if (v[j] >= id) break;
            }
        }
    }
    else if (flag == localTinyBitmap){

        bool rm_ok = true;

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

        TinyBitmap t_bmp(&setPtrTinyBmp);

        for (; (i < v.size()) && rm_ok; ++i) rm_ok = t_bmp.remove(v[i]);

        if (rm_ok){

            const size_t card = t_bmp.size();

            if (card == 0){

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                clear();
            }
            else if (card == 1){

                const uint32_t id = *(t_bmp.begin());

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;

                clear();
                add(id);
            }
            else if (t_bmp.maximum() < maxBitVectorIDs){

                PairID new_uc;

                for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) new_uc.add(*it);

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;

                *this = move(new_uc);
            }
            else setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        }
        else {

            const size_t sz_t_bmp = t_bmp.size();

            uint32_t* values = new uint32_t[sz_t_bmp];

            size_t j = 0;

            Bitmap* setPtrBmp = new Bitmap;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) values[j++] = *it;

            t_bmp.clear();
            setPtrBmp->r.addMany(sz_t_bmp, values);

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) {

        Bitmap* bitmap = getPtrBitmap();

        for (; i < v.size(); ++i) bitmap->r.remove(v[i]);

        const size_t card = bitmap->r.cardinality();

        if (card == 0) clear();
        else if (card == 1){

            const uint32_t id = bitmap->r.minimum();

            clear();
            add(id);
        }
        else if (bitmap->r.maximum() < maxBitVectorIDs){

            PairID new_uc;

            for (const_iterator it = begin(), it_end = end(); it != it_end; ++it) new_uc.add(*it);

            *this = move(new_uc);
        }
        else if ((setBits & flagMask) == ptrBitmap) bitmap->r.runOptimize();
    }
}

bool PairID::contains(const size_t pair_id) const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.contains(pair_id);

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const bool ret = t_bmp.contains(pair_id);

        t_bmp.detach();

        return ret;
    }

    if (flag == localSingleInt) return (pair_id == (setBits >> shiftMaskBits));

    if (pair_id < maxBitVectorIDs){

        const uintptr_t setBits_tmp = 0x1ULL << (pair_id + shiftMaskBits);

        return ((setBits & setBits_tmp) != 0);
    }

    return false;
}

size_t PairID::maximum() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.maximum();

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.maximum();

        t_bmp.detach();

        return ret;
    }

    if (flag == localSingleInt) return (setBits >> shiftMaskBits);

    const int nb_lead_0 = __builtin_clzll(setBits | flagMask);

    return (maxBitVectorIDs - nb_lead_0 - (nb_lead_0 != maxBitVectorIDs));
}

size_t PairID::minimum() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.minimum();

    const_iterator it = begin(), it_end = end();

    if (it != it_end) return *it;
    return 0;
}

size_t PairID::size() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.cardinality();

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.size();

        t_bmp.detach();

        return ret;
    }

    if (flag == localBitVector) return __builtin_popcountll(setBits & pointerMask);

    return 1;
}

bool PairID::write(ostream& stream_out) const {

    if (stream_out.good()){

        const uintptr_t flag = setBits & flagMask;

        if (flag == ptrBitmap){

            const uint32_t expected_sz = getConstPtrBitmap()->r.getSizeInBytes();

            const uintptr_t flag_expected_sz = (static_cast<uintptr_t>(expected_sz) << shiftMaskBits) | flag;

            char* serialized = new char[expected_sz];

            getConstPtrBitmap()->r.write(serialized);

            stream_out.write(reinterpret_cast<const char*>(&flag_expected_sz), sizeof(uintptr_t));
            stream_out.write(serialized, expected_sz);

            delete[] serialized;
        }
        else if (flag == localTinyBitmap){

            uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
            TinyBitmap t_bmp(&setPtrTinyBmp);

            stream_out.write(reinterpret_cast<const char*>(&flag), sizeof(uintptr_t));

            t_bmp.write(stream_out);
            t_bmp.detach();
        }
        else stream_out.write(reinterpret_cast<const char*>(&setBits), sizeof(uintptr_t));

        return true;
    }

    return false;
}

bool PairID::read(istream& stream_in) {

    if (stream_in.good()){

        clear();

        stream_in.read(reinterpret_cast<char*>(&setBits), sizeof(uintptr_t));

        const uintptr_t flag = setBits & flagMask;

        if (flag == ptrBitmap){

            Bitmap* setPtrBmp = new Bitmap;

            const uint32_t expected_sz = static_cast<uint32_t>(setBits >> shiftMaskBits);

            char* serialized = new char[expected_sz];

            stream_in.read(serialized, expected_sz);

            setPtrBmp->r = Roaring::read(serialized);

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;

            delete[] serialized;
        }
        else if (flag == localTinyBitmap){

            TinyBitmap t_bmp;

            t_bmp.read(stream_in);

            setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        }

        return true;
    }

    return false;
}

void PairID::runOptimize() {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) getPtrBitmap()->r.runOptimize();

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);

        t_bmp.runOptimize();

        setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
    }
}

PairID PairID::subsample(const uint64_t nb) const {

    const size_t card = cardinality();

    if (card <= nb) return *this;
    if (nb == 0) return PairID();

    const uintptr_t flag = setBits & flagMask;

    PairID pid;

    if (flag == ptrBitmap) {

        const Roaring r = getPtrBitmap()->r.subsample(nb);

        for (const auto id : r) pid.add(id);
    }
    else {

        PairID p_pos;

        std::random_device rd; // Seed
        std::default_random_engine generator(rd()); // Random number generator
        std::uniform_int_distribution<> distribution(0, card - 1); // Distribution on which to apply the generator

        for (size_t i = 0; i < nb; ++i) p_pos.add(distribution(generator));

        PairID::const_iterator it_pos = p_pos.begin(), it_pos_e = p_pos.end();
        PairID::const_iterator it = begin();

        for (size_t i = 0; (it_pos != it_pos_e); ++i, ++it) {

            if (i == *it_pos) {

                pid.add(*it);
                ++it_pos;
            }
        }
    }

    pid.runOptimize();

    return pid;
}

PairID::const_iterator PairID::begin() const {

    const_iterator it(this, true);

    return ++it;
}

PairID::const_iterator PairID::begin(const size_t id) const {

    const_iterator it(this, true);
    const_iterator it_e(this, false);

    if (isEmpty()) return it_e;

    if (it.flag == ptrBitmap) {

        const Roaring& roar = it.cs->getConstPtrBitmap()->r;

        it.it_roar.equalorlarger(id);

        if (it.it_roar != roar.end()) {

            it.ck_id = *(it.it_roar);
            it.it_setBits = roar.rank(it.ck_id) - 1;
        }
        else return it_e;
    }
    else {

        ++it;

        if ((it.flag == localTinyBitmap) && ((it.ck_id >> 16) != (id >> 16))) return it_e;

        while ((it != it_e) && (*it < id)) ++it;
    }

    return it;
}

PairID::const_iterator PairID::end() const {

    return const_iterator(this, false);
}

PairID::PairID_const_iterator::PairID_const_iterator() : it_roar(empty_roar.end()) {

    clear();
}

PairID::PairID_const_iterator::PairID_const_iterator( const PairID* cs_, const bool beg) :  cs(cs_), ck_id(0xffffffffffffffff), it_setBits(0xffffffffffffffff),
                                                                                            it_roar(empty_roar.end()) {

    flag = cs->setBits & flagMask;

    if (flag == ptrBitmap){

        it_roar = beg ? cs->getConstPtrBitmap()->r.begin() : cs->getConstPtrBitmap()->r.end();
        cs_sz = cs->size();
    }
    else if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = cs->getPtrTinyBitmap();

        t_bmp = &setPtrTinyBmp;
        it_t_bmp = beg ? t_bmp.begin() : t_bmp.end();
        cs_sz = cs->size();
    }
    else cs_sz = (flag == localSingleInt) ? 1 : maxBitVectorIDs;

    if (!beg) it_setBits = cs_sz;
}

PairID::PairID_const_iterator::PairID_const_iterator(const PairID_const_iterator& o) :  cs(o.cs), flag(o.flag), it_setBits(o.it_setBits), cs_sz(o.cs_sz),
                                                                                        ck_id(o.ck_id), it_roar(o.it_roar), it_t_bmp(o.it_t_bmp) {

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = cs->getPtrTinyBitmap();
        t_bmp = &setPtrTinyBmp;
    }
}

PairID::PairID_const_iterator::~PairID_const_iterator() {

    t_bmp.detach();
}

void PairID::PairID_const_iterator::clear() {

    cs = nullptr;
    flag = localBitVector;
    it_setBits = 0;
    it_roar = empty_roar.end();
    ck_id = 0xffffffffffffffff;
    cs_sz = 0;

    t_bmp.detach();
}

PairID::PairID_const_iterator& PairID::PairID_const_iterator::operator=(const PairID_const_iterator& o) {

    cs = o.cs;

    flag = o.flag;

    it_setBits = o.it_setBits;
    cs_sz = o.cs_sz;

    ck_id = o.ck_id;

    it_roar = o.it_roar;
    it_t_bmp = o.it_t_bmp;

    if (flag == localTinyBitmap){

        t_bmp.detach();

        uint16_t* setPtrTinyBmp = cs->getPtrTinyBitmap();
        t_bmp = &setPtrTinyBmp;
    }

    return *this;
}

PairID::PairID_const_iterator PairID::PairID_const_iterator::operator++(int) {

    PairID_const_iterator tmp(*this);
    operator++();
    return tmp;
}

PairID::PairID_const_iterator& PairID::PairID_const_iterator::operator++() {

    if (it_setBits != cs_sz){

        ++it_setBits;

        if (flag == ptrBitmap) {

            if (it_setBits != 0) ++it_roar;
            if (it_roar != cs->getConstPtrBitmap()->r.end()) ck_id = *it_roar;
        }
        else if (flag == localTinyBitmap) {

            if (it_setBits != 0) ++it_t_bmp;
            if (it_t_bmp != t_bmp.end()) ck_id = *it_t_bmp;
        }
        else if (flag == localBitVector){

            while (it_setBits < maxBitVectorIDs){

                if (((cs->setBits >> (it_setBits + shiftMaskBits)) & 0x1) != 0){

                    ck_id = it_setBits;
                    break;
                }

                ++it_setBits;
            }
        }
        else ck_id = cs->setBits >> shiftMaskBits;
    }

    return *this;
}

bool PairID::PairID_const_iterator::operator==(const PairID_const_iterator& o) const {

    if ((cs == o.cs) && (flag == o.flag) && (cs_sz == o.cs_sz)){

        if (flag == ptrBitmap) return (it_roar == o.it_roar);
        if (flag == localTinyBitmap) return (it_t_bmp == o.it_t_bmp);

        return (it_setBits == o.it_setBits);
    }

    return false;
}

bool PairID::PairID_const_iterator::operator!=(const PairID_const_iterator& o) const {

    return !operator==(o);
}

const size_t PairID::maxBitVectorIDs = 61; // 64 bits - 3 bits for the color set type

const uintptr_t PairID::localTinyBitmap = 0x0;
const uintptr_t PairID::localBitVector = 0x1;
const uintptr_t PairID::localSingleInt = 0x2;
const uintptr_t PairID::ptrBitmap = 0x3;

const size_t PairID::shiftMaskBits = 3;

const uintptr_t PairID::flagMask = 0x7;
const uintptr_t PairID::pointerMask = 0xFFFFFFFFFFFFFFF8;
