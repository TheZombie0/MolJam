import numpy as np


class ScoringHelpersMixin:
    def calculate_quality_score(self, error_rate, max_score=10, threshold_low=10, threshold_high=50):
        """
        Quality scoring function that demonstrates the "barrel effect" - 
        as quality decreases, penalties become increasingly severe.
        """
        # Handle edge cases: perfect and worst quality
        if error_rate <= 0:
            return max_score
        if error_rate >= 100:
            return 0.0

        # Normalize error rate and thresholds to [0, 1]
        e = error_rate / 100.0
        T1 = threshold_low / 100.0
        T2 = threshold_high / 100.0

        # Segment 1: Low error rates (linear penalty)
        if error_rate < threshold_low:
            score = max_score * (1 - e)

        # Segment 2: Medium error rates (quadratic penalty)
        elif error_rate < threshold_high:
            range_size = T2 - T1
            c = max_score * (1 - T1)
            b = -max_score
            a = (max_score/4 - c - b*range_size) / (range_size**2)
            delta = e - T1
            score = a*(delta**2) + b*delta + c

        # Segment 3: High error rates (exponential penalty decay)
        else:
            range_size = T2 - T1
            c_val = max_score * (1 - T1)
            b_val = -max_score
            a_val = (max_score/4 - c_val - b_val*range_size) / (range_size**2)
            derivative_at_T2 = 2 * a_val * range_size + b_val
            A = max_score / 4
            k = derivative_at_T2 / A
            delta = e - T2
            score = A * np.exp(k * delta)

        return max(0.0, min(max_score, score))
    def score_low_count_issues(self, issue_count, dataset_size, max_score=10, 
                               count_thresholds=[0, 1, 5, 10, 20, 50], 
                               rate_threshold=0.01):
        """
        Custom scoring function for issues with typically low counts.
        Considers both absolute count and relative rate.
        """
        if issue_count == 0:
            return max_score
        
        # Calculate rate
        rate = issue_count / dataset_size if dataset_size > 0 else 0
        
        # Score based on absolute count (60% weight)
        count_score = max_score
        for i, threshold in enumerate(count_thresholds):
            if issue_count <= threshold:
                # Linear interpolation between thresholds
                if i == 0:
                    count_score = max_score
                else:
                    prev_threshold = count_thresholds[i-1]
                    score_range = max_score / len(count_thresholds)
                    score_at_prev = max_score - (i-1) * score_range
                    score_at_curr = max_score - i * score_range
                    
                    if threshold > prev_threshold:
                        ratio = (issue_count - prev_threshold) / (threshold - prev_threshold)
                        count_score = score_at_prev - ratio * (score_at_prev - score_at_curr)
                break
        else:
            # Beyond highest threshold
            count_score = max_score / len(count_thresholds)
        
        # Score based on rate (40% weight)
        if rate <= rate_threshold:
            # Use logarithmic scale for small rates
            if rate > 0:
                log_rate = np.log10(rate * 1000)  # Scale up for better differentiation
                rate_score = max_score * (1 - max(0, min(1, (log_rate + 3) / 6)))
            else:
                rate_score = max_score
        else:
            # Linear penalty for higher rates
            rate_score = max_score * (1 - min(1, rate / 0.1))
        
        # Combine scores
        final_score = 0.6 * count_score + 0.4 * rate_score
        
        return max(0.0, min(max_score, final_score))
    def score_count_based_issues(self, issue_count, dataset_size, max_score=10,
                                 severity='medium'):
        """
        Scoring function for count-based issues like duplicates and contradictions.
        Uses different severity levels.
        """
        if issue_count == 0:
            return max_score
        
        # Calculate rate
        rate = issue_count / dataset_size if dataset_size > 0 else 0
        
        # Define severity parameters
        severity_params = {
            'low': {'base_penalty': 0.05, 'rate_multiplier': 20, 'max_penalty': 0.5},
            'medium': {'base_penalty': 0.1, 'rate_multiplier': 50, 'max_penalty': 0.7},
            'high': {'base_penalty': 0.2, 'rate_multiplier': 100, 'max_penalty': 0.9}
        }
        
        params = severity_params.get(severity, severity_params['medium'])
        
        # Calculate penalty using logarithmic scale for count
        if issue_count > 0:
            log_penalty = np.log10(issue_count + 1) * params['base_penalty']
        else:
            log_penalty = 0
        
        # Add rate-based penalty
        rate_penalty = rate * params['rate_multiplier'] * params['base_penalty']
        
        # Total penalty with maximum cap
        total_penalty = min(params['max_penalty'], log_penalty + rate_penalty)
        
        score = max_score * (1 - total_penalty)
        
        return max(0.0, min(max_score, score))

