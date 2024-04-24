#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List get_prob_list(int n_p_h, NumericVector p0, int n1, int n11, int n2, int n21) {
    List out(n_p_h);
    for (int i = 0; i < n_p_h; ++i) {
        double p_i = p0[i];
        NumericMatrix matr_i(n11, n21);
        for (int j = 0; j < n11; ++j) {
            for (int k = 0; k < n21; ++k) {
                matr_i(j, k) = R::dbinom(j, n1, p_i, 0) * R::dbinom(k, n2, p_i, 0);
            }
        }
        out[i] = matr_i;
    }
    return out;
}

int which_min(NumericVector x, int y) {
    int out = 0;
    double m0 = INFINITY;
    for (int i = 0; i < y; ++i) {
        double m = x[i];
        if (m < m0) {
            m0 = m;
            out = i;
        }
    }
    return out;
}

// [[Rcpp::export]]
int which_loc(int loc_i, int loc_j, IntegerVector select_i, IntegerVector select_j) {
    int l = select_i.length();
    int out = NA_INTEGER;
    for (int i = 0; i < l; ++i) {
        if (loc_i == select_i[i]) {
            if (loc_j == select_j[i]) {
                out = i + 1;
                break;
            }
        }
    }
    return out;
}

bool in_loc(int loc_i, int loc_j, std::vector<int> select_i, std::vector<int> select_j, int l) {
    bool out = false;
    for (int i = l - 1; i > -1; --i) {
        if (loc_i == select_i[i]) {
            if (loc_j == select_j[i]) {
                out = true;
                break;
            }
        }
    }
    return out;
}

void get_pmore(int* p_more, int* more_loc_i1, int* more_loc_j1, int* more_loc_i2, int* more_loc_j2,
    int select_i_new, int select_j_new, int n1, int n2,
    std::vector<int> select_i, std::vector<int> select_j, bool twoside, int k) {
    int select_i_new1 = select_i_new - 1;
    int select_j_new2 = select_j_new + 1;

    bool use1 = false;
    if (select_i_new1 >= 0) {
        if (twoside) {
            use1 = select_i_new1 * n2 > select_j_new * n1;
        }
        else {
            use1 = true;
        }
        if (use1) {
            if (select_j_new == 0) {
                use1 = true;
            }
            else if (in_loc(select_i_new1, select_j_new - 1, select_i, select_j, k)) {
                use1 = true;
            }
            else {
                use1 = false;
            }
        }
    }

    bool use2 = false;
    if (select_j_new2 <= n2) {
        if (twoside) {
            use2 = select_i_new * n2 > select_j_new2 * n1;
        }
        else {
            use2 = true;
        }
        if (use2) {
            if (select_i_new == n1) {
                use2 = true;
            }
            else if (in_loc(select_i_new + 1, select_j_new2, select_i, select_j, k)) {
                use2 = true;
            }
            else {
                use2 = false;
            }
        }
    }
    if (use1 && use2) {
        *p_more = 2;
        *more_loc_i1 = select_i_new1;
        *more_loc_j1 = select_j_new;
        *more_loc_i2 = select_i_new;
        *more_loc_j2 = select_j_new2;
    }
    else if (use1) {
        *p_more = 1;
        *more_loc_i1 = select_i_new1;
        *more_loc_j1 = select_j_new;
    }
    else if (use2) {
        *p_more = 1;
        *more_loc_i1 = select_i_new;
        *more_loc_j1 = select_j_new2;
    }
    return;
}

// [[Rcpp::export]]
NumericVector BNDtest_get_prob(int loc_i, int loc_j,
    int n1, int n2, int n_p, int n_pm1, int n_p_h, List prob_list, bool twoside) {
    int loc_i_rev = n1 - loc_i;
    int loc_j_rev = n2 - loc_j;
    if (twoside) {
        NumericVector p(n_p_h);
        for (int i = 0; i < n_p_h; ++i) {
            NumericMatrix p_i = prob_list[i];
            p[i] = p_i(loc_i, loc_j) + p_i(loc_i_rev, loc_j_rev);
        }
        return p;
    }
    else {
        NumericVector p(n_p);
        for (int i = 0; i < n_p_h; ++i) {
            NumericMatrix p_i = prob_list[i];
            p[i] = p_i(loc_i, loc_j);
            p[n_pm1 - i] = p_i(loc_i_rev, loc_j_rev);
        }
        return p;
    }
}

// [[Rcpp::export]]
NumericVector path_get_prob(IntegerVector select_i, IntegerVector select_j, int loc_des,
    int n1, int n2, int n_p, int n_pm1, int n_p_h, List prob_list, bool twoside, bool trace) {
    std::string e;
    NumericVector out = BNDtest_get_prob(select_i[0], select_j[0], n1, n2, n_p, n_pm1, n_p_h, prob_list, twoside);
    if (trace) {
        e = std::string(" / ") + std::to_string(loc_des) + std::string("\r");
        std::string prt = std::to_string(1) + e;
        Rcout << prt;
    }
    if (loc_des > 1) {
        for (int i = 1; i < loc_des; ++i) {
            out += BNDtest_get_prob(select_i[i], select_j[i], n1, n2, n_p, n_pm1, n_p_h, prob_list, twoside);
            if (trace) {
                std::string prt = std::to_string(i + 1) + e;
                Rcout << prt;
            }
        }
    }
    return out;
}

// [[Rcpp::export]]
List BNDtest_loop(NumericVector p_prob_max, int p_l, IntegerVector p_i, IntegerVector p_j,
    int k, int n0, bool trace, int select_i0, int select_j0, List p_prob, List p_prob_indiv, int n1, int n2,
    std::vector<int> select_i, std::vector<int> select_j, bool twoside,
    int n_p, int n_pm1, int n_p_h, List prob_list, NumericVector p0, std::string alternative) {
    double p_value = -1;
    std::string e = std::string(" / ") + std::to_string(n0) + std::string("\r");
    while (true) {
        int loc = which_min(p_prob_max, p_l);
        int select_i_new = p_i[loc];
        int select_j_new = p_j[loc];
        if (trace) {
            std::string prt = std::to_string(k + 1) + e;
            Rcout << prt;
        }
        if (select_i_new == select_i0 && select_j_new == select_j0) {
            p_value = p_prob_max[loc];
            break;
        }
        NumericVector prob_work = p_prob[loc];

        --p_l;
        p_prob_indiv[loc] = p_prob_indiv[p_l];
        p_prob[loc] = p_prob[p_l];
        p_i[loc] = p_i[p_l];
        p_j[loc] = p_j[p_l];
        p_prob_indiv[p_l] = R_NilValue;
        p_prob[p_l] = R_NilValue;
        p_i[p_l] = NA_INTEGER;
        p_j[p_l] = NA_INTEGER;
        p_prob_max[p_l] = NA_REAL;
        if (p_l > 0) {
            for (int i = 0; i < p_l; ++i) {
                NumericVector prob_i0 = p_prob_indiv[i];
                NumericVector prob_i = prob_work + prob_i0;
                p_prob_max[i] = max(prob_i);
                p_prob[i] = prob_i;
            }
        }
        int p_more = 0;
        int more_loc_i1, more_loc_j1;
        int more_loc_i2, more_loc_j2;
        get_pmore(&p_more, &more_loc_i1, &more_loc_j1, &more_loc_i2, &more_loc_j2,
            select_i_new, select_j_new, n1, n2, select_i, select_j, twoside, k);
        if (p_more > 0) {
            NumericVector prob_i1 = BNDtest_get_prob(more_loc_i1, more_loc_j1,
                n1, n2, n_p, n_pm1, n_p_h, prob_list, twoside);
            p_prob_indiv[p_l] = prob_i1;
            NumericVector prob_i2 = prob_work + prob_i1;
            p_prob_max[p_l] = max(prob_i2);
            p_prob[p_l] = prob_i2;
            p_i[p_l] = more_loc_i1;
            p_j[p_l] = more_loc_j1;
            ++p_l;
            if (p_more == 2) {
                NumericVector prob_i1 = BNDtest_get_prob(more_loc_i2, more_loc_j2,
                    n1, n2, n_p, n_pm1, n_p_h, prob_list, twoside);
                p_prob_indiv[p_l] = prob_i1;
                NumericVector prob_i2 = prob_work + prob_i1;
                p_prob_max[p_l] = max(prob_i2);
                p_prob[p_l] = prob_i2;
                p_i[p_l] = more_loc_i2;
                p_j[p_l] = more_loc_j2;
                ++p_l;
            }
        }
        select_i.push_back(select_i_new);
        select_j.push_back(select_j_new);
        ++k;
    }

    List data = List::create(
        Named("n1") = n1, Named("n2") = n2, Named("p0") = p0, Named("alternative") = alternative, Named("prob_list") = prob_list,
        Named("k") = k, Named("p_l") = p_l, Named("p_prob_indiv") = p_prob_indiv, Named("p_prob") = p_prob, 
        Named("p_prob_max") = p_prob_max, Named("p_i") = p_i, Named("p_j") = p_j, Named("select_i") = wrap(select_i), Named("select_j") = wrap(select_j)
    );
    
    List out = List::create(Named("p.value") = p_value,
        Named("data") = data);
    
    return out;
}

// [[Rcpp::export]]
List copy_list(List x, int n, int n_each) {
    int l = x.length();
    List out(l);
    for (int i = 0; i < n; ++i) {
        NumericVector out_i(n_each);
        NumericVector x_i = x[i];
        for (int j = 0; j < n_each; ++j) {
            out_i[j] = x_i[j];
        }
        out[i] = out_i;
    }
    return out;
}
