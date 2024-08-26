import axios from 'axios';
import { jwtDecode } from 'jwt-decode';

const api = axios.create({
  baseURL: process.env.REACT_APP_BASE_URL,
  headers: {
    'Content-Type': 'application/x-www-form-urlencoded'
}
});


api.interceptors.request.use(
  (config) => {
    const token = localStorage.getItem('access_token');
    if (token) {
      config.headers.Authorization = `Bearer ${token}`;
    } else {
      localStorage.setItem('username', null)
    }
    return config;
  },
  (error) => Promise.reject(error)
);


api.interceptors.response.use(
  (response) => response,
  async (error) => {
    const originalRequest = error.config;

    if (error.response.status === 401 && !originalRequest._retry) {
      originalRequest._retry = true;

      try {
        const refreshToken = localStorage.getItem('refresh_token');
        const response = await axios.post('/auth/refresh', {
          refresh_token: refreshToken
        });
        const { token } = response.data.access_token;

        localStorage.setItem('access_token', token);

        let decoded_token = jwtDecode(token)
        localStorage.setItem('username', decoded_token.sub)

        originalRequest.headers.Authorization = `Bearer ${token}`;
        return axios(originalRequest);
      } catch (error) {
        throw error
      }
    }
    return Promise.reject(error);
  }
);

export default api;
