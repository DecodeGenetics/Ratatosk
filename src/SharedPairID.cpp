#include "SharedPairID.hpp"

SharedPairID::SharedPairID() : global_pid(reinterpret_cast<uintptr_t>(nullptr)) {

	clear();
}

SharedPairID::SharedPairID(const SharedPairID& o) : local_pid(o.local_pid) {

	const PairID* g_pid = o.getGlobalPairID();

	if ((g_pid == nullptr) || g_pid->isEmpty()) global_pid = reinterpret_cast<uintptr_t>(nullptr);
	else {

		PairID* pid = new PairID;

		*pid = *g_pid;
		global_pid = reinterpret_cast<uintptr_t>(pid) | ptr_owner;
	}
}

SharedPairID::SharedPairID(SharedPairID&& o) : local_pid(move(o.local_pid)) {

	global_pid = o.global_pid;
	o.global_pid = reinterpret_cast<uintptr_t>(nullptr);

	o.clear();
}

SharedPairID::SharedPairID(PairID&& pid) : global_pid(reinterpret_cast<uintptr_t>(nullptr)), local_pid(move(pid)) {

	pid.clear();
}

SharedPairID::SharedPairID(const PairID& pid) : global_pid(reinterpret_cast<uintptr_t>(nullptr)), local_pid(pid) {}

SharedPairID::~SharedPairID() {

    clear();
}

void SharedPairID::clear() {

	if (hasGlobalOwnership()) delete getGlobalPairID();

	global_pid = reinterpret_cast<uintptr_t>(nullptr);

	local_pid.clear();
}

SharedPairID& SharedPairID::operator=(const SharedPairID& o){

    if (this != &o) {

    	const PairID* g_pid = o.getGlobalPairID();

    	clear();

    	local_pid = o.local_pid;

		if ((g_pid != nullptr) && !g_pid->isEmpty()) {

			PairID* pid = new PairID;

			*pid = *g_pid;
			global_pid = reinterpret_cast<uintptr_t>(pid) | ptr_owner;
		}
    }

    return *this;
}

SharedPairID& SharedPairID::operator=(SharedPairID&& o){

    if (this != &o) {

    	clear();

    	local_pid = move(o.local_pid);
		global_pid = o.global_pid;

		o.global_pid = reinterpret_cast<uintptr_t>(nullptr);

		o.clear();
    }

    return *this;
}

/*SharedPairID& SharedPairID::operator=(PairID&& pid){

    clear();

   	local_pid = move(pid);

   	pid.clear();

    return *this;
}

SharedPairID& SharedPairID::operator=(const PairID& pid){

	clear();

	local_pid = pid;

    return *this;
}*/

bool SharedPairID::operator==(const SharedPairID& o) const {

    if (size() != o.size()) return false;

    SharedPairID::const_iterator it(begin()), o_it(o.begin());
    const SharedPairID::const_iterator it_end(end()), o_it_end(o.end());

    for (; (it != it_end) && (o_it != o_it_end); ++it, ++o_it){

        if (*it != *o_it) return false;
    }

    return ((it == it_end) && (o_it == o_it_end));
}

SharedPairID SharedPairID::operator|(const PairID& rhs) const {

    SharedPairID lhs(*this);

    lhs |= rhs;
    return lhs;
}

SharedPairID SharedPairID::operator|(const SharedPairID& rhs) const {

    SharedPairID lhs(*this);

    lhs |= rhs;
    return lhs;
}

SharedPairID& SharedPairID::operator|=(const PairID& rhs) {

    if (cardinality() == 0) setLocalPairID(rhs);
    else if (rhs.cardinality() != 0){

		const PairID* g_pid = getGlobalPairID();

		if ((g_pid == nullptr) || g_pid->isEmpty()) local_pid |= rhs;
		else {

			const size_t a_card = cardinality();
			const size_t b_card = rhs.cardinality();

			const size_t log2_a = b_card * l_approximate_log2(a_card);
			const size_t log2_b = a_card * l_approximate_log2(b_card);

			const size_t min_a_b = min(a_card + b_card, min(log2_a, log2_b));

			PairID pid;
			
			if (min_a_b == log2_a) {

				PairID::const_iterator b_it_s = rhs.begin(), b_it_e = rhs.end();

				while (b_it_s != b_it_e){

					if (!contains(*b_it_s)) pid.add(*b_it_s);

					++b_it_s;
				}
			}
			else {

				SharedPairID::const_iterator a_it_s = begin(), a_it_e = end();
				PairID::const_iterator b_it_s = rhs.begin(), b_it_e = rhs.end();

				while ((a_it_s != a_it_e) && (b_it_s != b_it_e)){

					const uint32_t val_a = *a_it_s;
					const uint32_t val_b = *b_it_s;

					if (val_a == val_b){

						++a_it_s;
						++b_it_s;
					}
					else if (val_a < val_b) ++a_it_s;
					else {

						pid.add(val_b);

						++b_it_s;
					}
				}

				while (b_it_s != b_it_e) {

					pid.add(*b_it_s);

					++b_it_s;
				}
			}

			local_pid |= pid;
		}
    }

    return *this;
}

SharedPairID& SharedPairID::operator|=(const SharedPairID& rhs) {

    if (cardinality() == 0) *this = rhs;
    else if (rhs.cardinality() != 0){

		const PairID* gpid = getGlobalPairID();
		const PairID* rhs_gpid = rhs.getGlobalPairID();

		if ((gpid == nullptr) || gpid->isEmpty()) {

			local_pid |= rhs.local_pid;

			if ((rhs_gpid != nullptr) && !rhs_gpid->isEmpty()) local_pid |= *rhs_gpid;
		}
		else if (rhs_gpid == gpid) local_pid |= rhs.local_pid;
		else {

			const size_t a_card = cardinality();
			const size_t b_card = rhs.cardinality();

			const size_t log2_a = b_card * l_approximate_log2(a_card);
			const size_t log2_b = a_card * l_approximate_log2(b_card);

			const size_t min_a_b = min(a_card + b_card, min(log2_a, log2_b));

			PairID pid;
			
			if (min_a_b == log2_a) {

				SharedPairID::const_iterator b_it_s = rhs.begin(), b_it_e = rhs.end();

				while (b_it_s != b_it_e){

					if (!contains(*b_it_s)) pid.add(*b_it_s);

					++b_it_s;
				}
			}
			else {

				SharedPairID::const_iterator a_it_s = begin(), a_it_e = end();
				SharedPairID::const_iterator b_it_s = rhs.begin(), b_it_e = rhs.end();

				while ((a_it_s != a_it_e) && (b_it_s != b_it_e)){

					const uint32_t val_a = *a_it_s;
					const uint32_t val_b = *b_it_s;

					if (val_a == val_b){

						++a_it_s;
						++b_it_s;
					}
					else if (val_a < val_b) ++a_it_s;
					else {

						pid.add(val_b);

						++b_it_s;
					}
				}

				while (b_it_s != b_it_e) {

					pid.add(*b_it_s);

					++b_it_s;
				}
			}

			local_pid |= pid;
		}
    }

    return *this;
}

size_t SharedPairID::getSizeInBytes() const {

    return sizeof(SharedPairID) + local_pid.getSizeInBytes() + (hasGlobalOwnership() ? getGlobalPairID()->getSizeInBytes() : 0);
}

void SharedPairID::add(const size_t id) {

	const PairID* g_pid = getGlobalPairID();

	if ((g_pid == nullptr) || !g_pid->contains(id)) local_pid.add(id);
}

bool SharedPairID::contains(const size_t id) const {

	const PairID* g_pid = getGlobalPairID();

	return (local_pid.contains(id) || ((g_pid != nullptr) && g_pid->contains(id)));
}

size_t SharedPairID::maximum() const {

	const PairID* g_pid = getGlobalPairID();

	if (g_pid != nullptr) return max(local_pid.maximum(), g_pid->maximum());

	return local_pid.maximum();
}

size_t SharedPairID::minimum() const {

	const PairID* g_pid = getGlobalPairID();

	if (g_pid != nullptr) return min(local_pid.minimum(), g_pid->minimum());

	return local_pid.minimum();
}

size_t SharedPairID::size() const {

	const PairID* g_pid = getGlobalPairID();

    return local_pid.size() + ((g_pid != nullptr) ? g_pid->size() : 0);
}

void SharedPairID::runOptimize() {

    local_pid.runOptimize();
}

PairID SharedPairID::toPairID() const {

	if (getGlobalPairID() == nullptr) return local_pid;

	return (*getGlobalPairID() | local_pid);
}

vector<uint32_t> SharedPairID::toVector() const {

    vector<uint32_t> v;

    v.reserve(cardinality());

    for (const uint32_t id : *this) v.push_back(id);

    return v;
}

pair<const PairID*, const PairID*> SharedPairID::getPairIDs() const {

	return pair<const PairID*, const PairID*>(getGlobalPairID(), &local_pid);
}

const PairID& SharedPairID::getLocalPairID() const {

	return local_pid;
}

void SharedPairID::setGlobalPairID(const PairID* g_pid) {

	const PairID l_pid = move(local_pid);

	clear();

	global_pid = reinterpret_cast<uintptr_t>(g_pid) | ptr_user;

	setLocalPairID(l_pid);
}

void SharedPairID::setGlobalPairID(const PairID& g_pid) {

	const PairID l_pid = move(local_pid);

	clear();

	PairID* pid = new PairID;

	*pid = g_pid;
	global_pid = reinterpret_cast<uintptr_t>(pid) | ptr_owner;

	setLocalPairID(l_pid);
}

void SharedPairID::setLocalPairID(const PairID& l_pid) {

	const PairID* g_pid = getGlobalPairID();

	if ((g_pid == nullptr) || g_pid->isEmpty() || l_pid.isEmpty()) local_pid = l_pid;
	else {

		const size_t a_card = g_pid->cardinality();
		const size_t b_card = l_pid.cardinality();

		const size_t log2_a = b_card * l_approximate_log2(a_card);
		const size_t log2_b = a_card * l_approximate_log2(b_card);

		const size_t min_a_b = min(a_card + b_card, min(log2_a, log2_b));

		PairID pid;

		local_pid.clear();
		
		if (min_a_b == log2_a) {

			PairID::const_iterator b_it_s = l_pid.begin(), b_it_e = l_pid.end();

			while (b_it_s != b_it_e){

				if (!g_pid->contains(*b_it_s)) pid.add(*b_it_s);

				++b_it_s;
			}
		}
		else {

			PairID::const_iterator a_it_s = g_pid->begin(), a_it_e = g_pid->end();
			PairID::const_iterator b_it_s = l_pid.begin(), b_it_e = l_pid.end();

			while ((a_it_s != a_it_e) && (b_it_s != b_it_e)){

				const uint32_t val_a = *a_it_s;
				const uint32_t val_b = *b_it_s;

				if (val_a == val_b){

					pid.add(val_b);

					++a_it_s;
					++b_it_s;
				}
				else if (val_a < val_b) ++a_it_s;
				else {

					pid.add(val_b);

					++b_it_s;
				}
			}

			while (b_it_s != b_it_e) {

				pid.add(*b_it_s);

				++b_it_s;
			}
		}

		local_pid = move(pid);
	}
}

bool SharedPairID::write(ostream& stream_out) const {

	const PairID* g_pid = getGlobalPairID();

	if (g_pid != nullptr) {

		if (stream_out.good()) g_pid->write(stream_out);
	}
	else {

		PairID dummy_pid;

		if (stream_out.good()) dummy_pid.write(stream_out);
	}
	
	if (stream_out.good()) local_pid.write(stream_out);

	return (stream_out.good());
}

bool SharedPairID::read(istream& stream_in) {

	PairID g_pid;

    clear();

    g_pid.read(stream_in); // Read global PairID

    if (!g_pid.isEmpty()) setGlobalPairID(g_pid);

    local_pid.read(stream_in);

	return (stream_in.good());
}

SharedPairID::const_iterator SharedPairID::begin() const {

    return const_iterator(this, true);
}

SharedPairID::const_iterator SharedPairID::end() const {

    return const_iterator(this, false);
}

SharedPairID::SharedPairID_const_iterator::SharedPairID_const_iterator() {

	clear();
}

SharedPairID::SharedPairID_const_iterator::SharedPairID_const_iterator(const SharedPairID* spid, const bool begin) {

	if (spid == nullptr) clear();
	else {

		{
			local_pid = &(spid->local_pid);

			if (begin) local_pid_it = local_pid->begin();
			else local_pid_it = local_pid->end();

			l_valid = (local_pid_it != local_pid->end());
		}

		{
			global_pid = spid->getGlobalPairID();

			if (global_pid == nullptr) {

				global_pid = local_pid;
				global_pid_it = global_pid->end();
				g_valid = false;
			}
			else {

				if (begin) global_pid_it = global_pid->begin();
				else global_pid_it = global_pid->end();

				g_valid = (global_pid_it != global_pid->end());
			}
		}
	}
}

SharedPairID::SharedPairID_const_iterator::SharedPairID_const_iterator(const SharedPairID_const_iterator& o) : 	global_pid(o.global_pid), local_pid(o.local_pid),
																												global_pid_it(o.global_pid_it), local_pid_it(o.local_pid_it),
																												g_valid(o.g_valid), l_valid(o.l_valid) {}

void SharedPairID::SharedPairID_const_iterator::clear() {

	global_pid = nullptr;
	local_pid = nullptr;

	g_valid = false;
	l_valid = false;

	local_pid_it.clear();
	global_pid_it.clear();
}

uint32_t SharedPairID::SharedPairID_const_iterator::operator*() const {

    if (l_valid && g_valid) return min(*local_pid_it, *global_pid_it);
    else if (g_valid) return *global_pid_it;
    
    return *local_pid_it;
}

SharedPairID::SharedPairID_const_iterator& SharedPairID::SharedPairID_const_iterator::operator=(const SharedPairID_const_iterator& o) {

	if (this != &o) {

		global_pid = o.global_pid;
		local_pid = o.local_pid;

		global_pid_it = o.global_pid_it;
		local_pid_it = o.local_pid_it;

		g_valid = o.g_valid;
		l_valid = o.l_valid;
	}

    return *this;
}

SharedPairID::SharedPairID_const_iterator SharedPairID::SharedPairID_const_iterator::operator++(int) {

    SharedPairID_const_iterator tmp(*this);

    operator++();

    return tmp;
}

SharedPairID::SharedPairID_const_iterator& SharedPairID::SharedPairID_const_iterator::operator++() {

	if (l_valid || g_valid) {

		if (l_valid && g_valid) {

			if (*local_pid_it < *global_pid_it) ++local_pid_it;
			else if (*local_pid_it > *global_pid_it) ++global_pid_it;
			else {

				++local_pid_it;
				++global_pid_it;
			}
		}
		else if (l_valid) ++local_pid_it;
		else ++global_pid_it;

		l_valid = (local_pid_it != local_pid->end());
		g_valid = (global_pid_it != global_pid->end());
	}

    return *this;
}

bool SharedPairID::SharedPairID_const_iterator::operator==(const SharedPairID_const_iterator& o) const {

    return (local_pid_it == o.local_pid_it) && (global_pid_it == o.global_pid_it) && (l_valid == o.l_valid) && (g_valid == o.g_valid);
}

bool SharedPairID::SharedPairID_const_iterator::operator!=(const SharedPairID_const_iterator& o) const {

    return !operator==(o);
}

const uintptr_t SharedPairID::ptr_user = reinterpret_cast<uintptr_t>(nullptr) & 0x1; // This implies that we are never owner of a nullptr
const uintptr_t SharedPairID::ptr_owner = ~reinterpret_cast<uintptr_t>(nullptr) & 0x1;

const size_t SharedPairID::shiftMaskBits = 1;

const uintptr_t SharedPairID::flagMask = 0x1ULL;
const uintptr_t SharedPairID::ptrMask = 0xfffffffffffffffeULL;