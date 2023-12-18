#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  ResultVal status = EXIT_OK;
  if (rows < 1 || columns < 1 || !result)
    status = EXIT_INCORR_MATRIX;
  else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (!result->matrix) {
      s21_remove_matrix(result);
      status = EXIT_INCORR_MATRIX;
    } else {
      for (int i = 0; i < rows; ++i) {
        result->matrix[i] = (double *)calloc(columns, sizeof(double));
        if (!result->matrix[i]) {
          s21_remove_matrix(result);
          status = EXIT_INCORR_MATRIX;
          break;
        }
      }
    }
  }

  return status;
}

void s21_remove_matrix(matrix_t *A) {
  if (A && A->matrix) {
    for (int i = 0; i < A->rows; ++i) {
      if (A->matrix[i]) free(A->matrix[i]);
    }
    free(A->matrix);
    A->matrix = NULL;
  }
  if (A && A->columns) A->columns = 0;
  if (A && A->rows) A->rows = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = SUCCESS;
  int equality_status = 1;  // initially, the elements of the matrix are equal
  if (s21_is_not_correct_matrix(A) || s21_is_not_correct_matrix(B))
    status = FAILURE;
  else {
    if (!s21_has_equal_sizes(A, B))
      status = FAILURE;
    else {
      for (int i = 0; i < A->rows; ++i) {
        for (int j = 0; j < A->columns; ++j) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1.0E-7) {
            status = FAILURE;
            equality_status = 0;
            break;
          }
        }
        if (!equality_status) break;
      }
    }
  }

  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  ResultVal status = EXIT_OK;
  if (s21_is_not_correct_matrix(A) || s21_is_not_correct_matrix(B))
    status = EXIT_INCORR_MATRIX;
  else if (!s21_has_equal_sizes(A, B))
    status = EXIT_CALC_ERROR;
  else {
    status = s21_create_matrix(A->rows, A->columns, result);
    if (!status) {
      for (int i = 0; i < A->rows; ++i)
        for (int j = 0; j < A->columns; ++j)
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }

  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  ResultVal status = EXIT_OK;
  if (s21_is_not_correct_matrix(A) || s21_is_not_correct_matrix(B))
    status = EXIT_INCORR_MATRIX;
  else if (!s21_has_equal_sizes(A, B))
    status = EXIT_CALC_ERROR;
  else {
    status = s21_create_matrix(A->rows, A->columns, result);
    if (!status) {
      for (int i = 0; i < A->rows; ++i)
        for (int j = 0; j < A->columns; ++j)
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }

  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  ResultVal status = EXIT_OK;
  if (s21_is_not_correct_matrix(A))
    status = EXIT_INCORR_MATRIX;
  else {
    status = s21_create_matrix(A->rows, A->columns, result);
    if (!status) {
      for (int i = 0; i < A->rows; ++i)
        for (int j = 0; j < A->columns; ++j)
          result->matrix[i][j] = number * A->matrix[i][j];
    }
  }

  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  ResultVal status = EXIT_OK;
  if (s21_is_not_correct_matrix(A) || s21_is_not_correct_matrix(B))
    status = EXIT_INCORR_MATRIX;
  else if (A->columns != B->rows)
    status = EXIT_CALC_ERROR;
  else {
    status = s21_create_matrix(A->rows, B->columns, result);
    if (!status) {
      for (int i = 0; i < A->rows; ++i)
        for (int j = 0; j < B->columns; ++j)
          for (int k = 0; k < B->rows; ++k)
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
    }
  }

  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  ResultVal status = EXIT_OK;
  if (s21_is_not_correct_matrix(A))
    status = EXIT_INCORR_MATRIX;
  else {
    status = s21_create_matrix(A->columns, A->rows, result);
    if (!status) {
      for (int i = 0; i < A->rows; ++i)
        for (int j = 0; j < A->columns; ++j)
          result->matrix[j][i] = A->matrix[i][j];
    }
  }

  return status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  ResultVal status = EXIT_OK;
  if (A->rows != A->columns && A->matrix) {
    status = EXIT_CALC_ERROR;
  } else if (!s21_is_not_correct_matrix(A) &&
             !s21_create_matrix(A->rows, A->columns, result)) {
    s21_calc_minor(A, result);
    for (int i = 0; i < result->rows; ++i)
      for (int j = 0; j < result->columns; ++j)
        result->matrix[i][j] *= (i + j) & 1 ? -1 : 1;
    status = EXIT_OK;
  } else {
    status = EXIT_INCORR_MATRIX;
  }

  return status;
}

int s21_determinant(matrix_t *A, double *result) {
  ResultVal status = EXIT_OK;
  if (A->rows != A->columns && A->matrix)
    status = EXIT_CALC_ERROR;
  else if (!s21_is_not_correct_matrix(A)) {
    *result = s21_calc_det(A);
    status = EXIT_OK;
  } else
    status = EXIT_INCORR_MATRIX;

  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  ResultVal status = EXIT_OK;
  if (s21_is_not_correct_matrix(A)) {
    status = EXIT_INCORR_MATRIX;
  } else {
    double det_A;
    if (s21_determinant(A, &det_A) || det_A == 0) {
      status = EXIT_CALC_ERROR;
    } else {
      matrix_t M;
      matrix_t M_transpose;
      if (s21_calc_complements(A, &M) || s21_transpose(&M, &M_transpose)) {
        status = EXIT_CALC_ERROR;
      } else {
        if (s21_mult_number(&M_transpose, 1.0 / det_A, result)) {
          status = EXIT_CALC_ERROR;
        }
      }
      s21_remove_matrix(&M_transpose);
      s21_remove_matrix(&M);
    }
  }

  return status;
}